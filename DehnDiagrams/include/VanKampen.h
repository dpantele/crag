/**
 * Simple inmplementation of Dehn Diagrams - planar graphs representing identity in the FP group
 */
#pragma once
#ifndef CRAG_DEHNDIAGRAM_H
#define CRAG_DEHNDIAGRAM_H

#include <assert.h>
#include <deque>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <tuple>
#include <vector>

template<typename WordT>
class VanKampen {
 public:
  using Word = WordT;
  using Letter = typename Word::Letter;

  VanKampen()
   : vertices_(1)
  { }

  struct Vertex;

  struct Edge {
    Vertex* endpoint_;
    Letter label_;
    Edge* inverse_ = nullptr;
    Edge* next_;
    Edge* prev_;
    size_t id_;
    static size_t GetNextId() {
      static size_t id = 0;
      return ++id;
    }

    Edge(Vertex* endpoint, Letter label)
      : endpoint_(endpoint)
      , label_(label)
      , next_(this)
      , prev_(this)
      , id_(GetNextId())
    { }

    Edge(const Edge&) = delete;
    Edge(Edge&&) = delete;

    Edge* next() {
      return next_;
    }

    const Edge* next() const {
      return next_;
    }

    Edge* prev() {
      return prev_;
    }

    void set_next(Edge* next) {
      next_ = next;
      next->prev_ = this;
    }

    void Insert(Edge* next) {
      auto current_next = next_;
      set_next(next);
      next->set_next(current_next);
    }

    void RemoveSelf() {
      prev_->set_next(next_);
      next_ = prev_ = this;
    }
  };

  struct Vertex {
    size_t id_;
    static size_t GetNextId() {
      static size_t id = 0;
      return ++id;
    }

    Vertex()
      : id_(GetNextId())
    { }

    Edge* first_ = nullptr;

    void Insert(Edge* edge) {
      if (first_ == nullptr) {
        first_ = edge;
      } else {
        first_->prev()->Insert(edge);
      }
    }

    bool IsEmpty() const {
      return first_ == nullptr;
    }

    Edge* first() {
      return first_;
    }

    const Edge* first() const {
      return first_;
    }

    void Move(Vertex* other, Edge* this_first, Edge* other_first, Edge* other_last) {
      assert(other->first_); //other can't be empty
      //first, we remove other_first .. other_last from the 'other' edges list
      if (other_last->next() != other_first) {
        //if they do not take the whole set of edges
        other_first->prev()->set_next(other_last->next());
        //also always reassign other->first_ to other_last->next()
        other->first_ = other_last->next();
      } else {
        //if they are all edges in the other vertex, just make the other vertex empty
        other->first_ = nullptr;
      }

      //if user have not passed this_first parameter, then this vertex must be empty
      if (!this_first) {
        //in this case simply make a list from other_first .. other_last
        //and associate it with this
        assert(this->IsEmpty());
        this->first_ = other_first;
        this_first = other_first;
        other_last->set_next(other_first);
      } else {
        //if current vertex is not empty, glue other_first .. other_last before this_first
        this_first->prev()->set_next(other_first);
        other_last->set_next(this_first);
      }

      //finally, change endpoint_ of just-moved-vertices
      auto current_other = other_first;
      do {
        assert(current_other->inverse_->endpoint_ == other);
        current_other->inverse_->endpoint_ = this;
        //if it is a loop, then the other end should also be included in the range,
        //so we don't need the code below
//        if (current_other->endpoint_ == other) {
//          current_other->endpoint_ = this;
//        }
        current_other = current_other->next();
      } while (current_other != this_first);
    }

    //! Move edges from other to this (they will be 'prepended' so that
    // other_last will be followed by this_first)
    void Join(Vertex* other, Edge* this_first, Edge* other_last) {
      assert(this != other);
      assert(other->first_);
      assert(first_);

      Move(other, this_first, other_last->next(), other_last);
    }

    //! Private function which removes a single edge and returns true if there is something else in the vertex
    bool RemoveBasic(Edge* e) {
      if (e == first_) {
        if (e->next() == e) {
          first_ = nullptr;
          return false;
        } else {
          first_ = first_->next();
        }
      }
      //e != first_
      e->RemoveSelf();
      return true;
    }

    //! Remove edge e and return the next one
    Edge* Remove(Edge* e) {
      e->endpoint_->RemoveBasic(e->inverse_);
      auto next = e->next();
      if (RemoveBasic(e)) {
        return next;
      }
      return nullptr;
    }
  };

  Vertex* root() {
    return &vertices_.front();
  }

  void Push(Word w) {
    w.Reduce();

    if (w.IsTrivial()) {
      return;
    }
    auto current_vertex = root();

    //stem
    while (w.length() > 1 && w.front() == -w.back()) {
      current_vertex = PushEdge(current_vertex, w.front());
      w.PopFront();
      w.PopBack();
    }

    //lollipop, except the last edge
    auto lollipop_start = current_vertex;

    while (w.length() > 1) {
      current_vertex = PushEdge(current_vertex, w.front());
      w.PopFront();
    }

    //the last edge
    //it is ok to insert it like that because of the lollipop structure
    PushEdge(current_vertex, lollipop_start, w.front());
  }

  void Fold() {
    //initially we have only root to check
    std::set<Vertex*> to_check = {root()};

    while (!to_check.empty()) {
      //for description of the folding process, see chapter
      //"Van Kampen Diagrams and Pictures"
      //from "The Geometry of the Word Problem for Finitely Generated Groups"
      //"Advanced Courses in Mathematics CRM Barcelona", p. 163-177

      Vertex* v = *to_check.begin();
      to_check.erase(to_check.begin());

      if (v->IsEmpty()) {
        continue;
      }

      auto e1 = v->first();
      auto e2 = e1->next();

      //if there are at least two edges, they could possibly be folded
      while (e1 != e2) {
        if (e1->label_ == e2->label_) {
          //fold here

          auto v1 = e1->endpoint_;
          auto v2 = e2->endpoint_;
          //case 1:
          if (v1 != v2 && v2 != v && v != v1) {
            //simply merge v2 and v1 and merge e1 and e2
            v1->Join(
              v2,
              e1->inverse_,
              e2->inverse_
            );
            //remove e2
            e2 = v->Remove(e2);

            to_check.emplace(v1);
          } else if (v1 == v && v2 != v) {
            //e1 is a loop
            //and e2 is not
            //fold e2 onto the e1
            v1->Join(
              v2,
              e1->inverse_,
              e2->inverse_
            );
            //remove e2
            v->Remove(e2);

            //and we will have to check v again
            to_check.emplace(v1);
            //so break now
            break;
          } else if (v2 == v && v1 != v) {
            //e2 is a loop
            //and e1 is not
            v1->Join(
              v2,
              e1->inverse_,
              e2->inverse_
            );
            //remove e2
            v->Remove(e2);
            to_check.emplace(v1);
            //v is empty now, just break
            break;
          } else if (v1 == v2 && v != v1) {
            //e1 and e2 are parallel
            //but are not loops

            //make two diagrams - one for the interior and another for the outside

            //if interior is not empty
            if (e2->inverse_->next() != e1->inverse_) {
              //first, create a new vertex which would serve as a base for a new sphere
              auto S = NewVertex();

              //now move all edges from the interior to a new vertex
              S->Move(v1, nullptr, e2->inverse_->next(), e1->inverse_->prev());
              to_check.emplace(S);
            }

            assert(e2->inverse_->next() == e1->inverse_);

            //and now remove e2 itself
            e2 = v->Remove(e2);
            to_check.emplace(v1);
          } else if (v1 == v2 && v == v1) {
            //two loops, both ending up in the same vertex
            //we will glue a new sphere

            //to glue two loops, we just need to make sure that there is nothing in between
            if (e2->inverse_->next() != e1->inverse_) {
              //if there is something, just move that to a separate sphere
              auto S = NewVertex();
              S->Move(v1, nullptr, e2->inverse_->next(), e1->inverse_->prev());
              //note that e1 and e2 could not have been affected
              //although v->first_ is e1->inverse_ now
              to_check.emplace(S);
            }

            assert(e2->inverse_->next() == e1->inverse_);

            //and now remove e2 itself
            v->Remove(e2);

            //since v was changed, just reconsider it
            to_check.emplace(v);
            break;
          }
        } else {
          e1 = e2;
          e2 = e1->next();

          //if we went through all edges, break here
          if (e1 == v->first()) {
            break;
          }
        }
      }

      if (v->IsEmpty()) {
        continue;
      }

      //if vertex has degree 1, remove the remaining edge
      if (v->first() == v->first()->next()) {
        //if a vertex has a single edge,
        //then we could simply remove this spike

        to_check.emplace(v->first()->endpoint_);
        v->Remove(v->first());
      }
    }
  }

  void PrintTo(std::ostream* out) const {
    std::map<const Vertex*, uint> vertex_ids;

    for (auto&& v : vertices_) {
      vertex_ids.emplace(&v, vertex_ids.size());
    }

    for (auto&& v : vertices_) {
      (*out) << "v" << /*vertex_ids[&v]*/v.id_ << ": ";
      auto edge = v.first();
      if (!edge) {
        *out << "empty\n";
        continue;
      }

      do {
        std::cout << edge->id_ << ": " << vertex_ids[edge->endpoint_] << " " << edge->label_ << ", ";
        edge = edge->next();
      } while (edge != v.first());

      *out << "\n";
    }
  }
 private:
  std::deque<Vertex> vertices_;
  std::deque<Edge> edges_;

  Vertex* NewVertex() {
    vertices_.emplace_back();
    return &vertices_.back();
  }

  Edge* NewEdge(Vertex* endpoint, Letter l) {
    edges_.emplace_back(endpoint, l);
    return &edges_.back();
  }

  void PushEdge(Vertex* origin, Vertex* endpoint, Letter l) {
    auto edge = NewEdge(endpoint, l);
    auto inverse = NewEdge(origin, -l);
    edge->inverse_ = inverse;
    inverse->inverse_ = edge;
    origin->Insert(edge);
    endpoint->Insert(inverse);
  }

  Vertex* PushEdge(Vertex* origin, Letter l) {
    auto new_vertex = NewVertex();
    PushEdge(origin, new_vertex, l);
    return new_vertex;
  }
};


#endif //CRAG_DEHNDIAGRAM_H
