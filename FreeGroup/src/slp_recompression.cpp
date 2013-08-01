#include "slp_recompression.h"
#include <iterator>

namespace crag {
namespace slp {
namespace recompression {

Rule::Rule(std::initializer_list<RuleLetter> letters)
  : debug_id(0)
{
  for (auto& letter : letters) {
    letters_.push_back(letter);

    if (letter.is_nonterminal()) {
      letters_.back().nonterminal_rule_->register_inclusion(
          this,
          std::prev(letters_.end())
      );
    }
  }

  if (!letters_.empty()) {
    first_terminal_letter_ = letters_.front().first_terminal_ptr();
    last_terminal_letter_ = letters_.back().last_terminal_ptr();
  }

  assert(first_terminal_letter_ && last_terminal_letter_);
}

Rule::iterator Rule::pop_first_from_letter(Rule::iterator letter_position) {
  assert(!empty());
  assert(!letter_position->is_empty_nonterminal());

  if (!letter_position->is_nonterminal()) {
    return letter_position;
  }

  Rule* letter_rule = letter_position->nonterminal_rule();

  iterator popping_letter = letter_rule->begin();
  assert(!popping_letter->is_empty_nonterminal());

  if (popping_letter->is_nonterminal()) {
    popping_letter = letter_rule->pop_first_from_letter(popping_letter);
  }

  assert(!popping_letter->is_nonterminal());
  assert(popping_letter == letter_rule->begin());

  const RuleLetter& popped_letter = *(letter_rule->delete_letter(popping_letter));

  assert(letter_rule->empty() || !letter_rule->begin()->is_empty_nonterminal());

  if (!letter_rule->empty()) {
    letter_rule->first_terminal_letter_ = letter_rule->begin()->first_terminal_ptr();
  } else {
    letter_rule->first_terminal_letter_ = nullptr;
    letter_rule->first_terminal_letter_ = nullptr;
  }

  iterator inserted_letter;

  for (auto& occurence : letter_rule->nonterminal_index_) {
    assert(!occurence.rule_->empty());

    occurence.rule_->insert_popped_letter_left(occurence.letter_, popped_letter);

    iterator popped_letter_inserted = std::prev(occurence.letter_);

    if (letter_rule->empty()) {
      assert(occurence.letter_->is_empty_nonterminal());
      popped_letter_inserted = occurence.rule_->remove_empty_letter(occurence.letter_).first;
    }

    assert(occurence.letter_ != begin());
    assert(!popped_letter_inserted->is_nonterminal());

    if (occurence.rule_ == this && occurence.letter_ == letter_position) {
      inserted_letter = popped_letter_inserted;
    } else {
      assert(!popped_letter_inserted->is_empty_terminal());
    }
  }

  return inserted_letter;
}

Rule::iterator Rule::pop_last_from_letter(Rule::iterator letter_position) {
  assert(!empty());
  assert(!letter_position->is_empty_nonterminal());

  if (!letter_position->is_nonterminal()) {
    return letter_position;
  }

  Rule* letter_rule = letter_position->nonterminal_rule();

  iterator popping_letter = std::prev(letter_rule->end());
  assert(!popping_letter->is_empty_nonterminal());

  if (popping_letter->is_nonterminal()) {
    popping_letter = letter_rule->pop_last_from_letter(popping_letter);
  }

  assert(!popping_letter->is_nonterminal());
  assert(popping_letter == std::prev(letter_rule->end()));

  const RuleLetter& popped_letter = *(letter_rule->delete_letter(popping_letter));

  assert(letter_rule->empty() || !std::prev(letter_rule->end())->is_empty_nonterminal());

  if (!letter_rule->empty()) {
    letter_rule->last_terminal_letter_ = std::prev(letter_rule->end())->last_terminal_ptr();
  } else {
    letter_rule->first_terminal_letter_ = nullptr;
    letter_rule->first_terminal_letter_ = nullptr;
  }

  iterator inserted_letter;

  for (auto& occurence : letter_rule->nonterminal_index_) {
    assert(!occurence.rule_->empty());

    occurence.rule_->insert_popped_letter_right(occurence.letter_, popped_letter);

    iterator popped_letter_inserted = std::next(occurence.letter_);

    if (letter_rule->empty()) {
      assert(occurence.letter_->is_empty_nonterminal());
      popped_letter_inserted = occurence.rule_->remove_empty_letter(occurence.letter_).second;
    }

    assert(popped_letter_inserted != end());
    assert(!popped_letter_inserted->is_nonterminal());

    if (occurence.rule_ == this && occurence.letter_ == letter_position) {
      inserted_letter = popped_letter_inserted;
    } else {
      assert(!popped_letter_inserted->is_empty_terminal());
    }
  }

  return inserted_letter;
}

std::pair<Rule::iterator, Rule::iterator> Rule::remove_empty_letter(Rule::iterator position) {
  assert(!empty());
  assert(position->is_empty_nonterminal() || position->is_empty_terminal());

  auto position_after = std::next(position);
  auto position_before = std::prev(position);

  delete_letter(position);

  if (empty()) {
    first_terminal_letter_ = nullptr;
    last_terminal_letter_ = nullptr;
    return std::make_pair(end(), end());
  }

  if (position_after == letters_.begin()) {
    return std::make_pair(end(), std::move(position_after));
  } else if (position_before == std::prev(end())) {
    return std::make_pair(std::move(position_before), end());
  }

  return std::make_pair(std::move(position_before), std::move(position_after));
}

Rule::iterator Rule::compress_power(Rule::iterator position, TerminalId new_terminal) {
  assert(!position->is_nonterminal());
  assert(!position->is_empty_terminal());

  assert(!begin()->is_empty_nonterminal());
  assert(!std::prev(end())->is_empty_nonterminal());

  position->terminal_.id = new_terminal;
  position->terminal_.power = position->terminal_sgn();

  assert(position == begin() || !std::prev(position)->is_power_of(new_terminal));
  assert(std::next(position) == end() || !std::next(position)->is_power_of(new_terminal));

  return position;
}

Rule::iterator Rule::compress_pair(
    Rule::iterator first,
    Rule::iterator second,
    TerminalId new_terminal
) {
  assert(second == std::next(first));
  assert(!first->is_nonterminal() && !first->is_power() && !first->is_empty_nonterminal());
  assert(!second->is_nonterminal() && !second->is_power() && !second->is_empty_nonterminal());

  assert(!begin()->is_empty_nonterminal());

  assert(!std::prev(end())->is_empty_nonterminal());

  first->terminal_.id = abs(new_terminal);
  first->terminal_.power = new_terminal > 0 ? 1 : -1;
  delete_letter(second);

  return first;
}

Rule::iterator Rule::merge_letters(Rule::iterator current, Rule::iterator next) {
  assert(current->is_valid());
  assert(next->is_valid());
  assert(std::next(current) == next);
  assert(!current->is_nonterminal());
  assert(!current->is_empty_terminal());
  assert(!next->is_nonterminal());
  assert(!next->is_empty_terminal());

  assert(current->is_power_of(next->first_terminal_letter_id()));
  assert(next->is_power_of(current->last_terminal_letter_id()));

  current->terminal_.power += next->terminal_power();
  delete_letter(next);

  return current;
}

void Rule::insert_popped_letter_right(
    Rule::iterator letter_position,
    const RuleLetter& popped_letter)
{
  assert(!popped_letter.is_nonterminal());
  assert(!popped_letter.is_valid());
  assert(letter_position->is_nonterminal());

  auto position_after = std::next(letter_position);

  iterator inserted = letters_.emplace(position_after, popped_letter);
}

void Rule::insert_popped_letter_left(
    Rule::iterator letter_position,
    const RuleLetter& popped_letter)
{
  assert(!popped_letter.is_nonterminal());
  assert(!popped_letter.is_valid());
  assert(letter_position->is_nonterminal());

  iterator inserted = letters_.emplace(letter_position, popped_letter);
}

void JezRules::initialize(const Vertex& slp)
{
  auto acceptor = [this] (const inspector::InspectorTask& task) {
    return this->vertex_rules_.count(task.vertex) == 0;
    //true only if vertex is not visited yet
  };
  Inspector<inspector::Postorder, decltype(acceptor)> inspector(slp, acceptor);

  while (!inspector.stopped()) {
    if (inspector.vertex().height() < 2) {
      get_letter(inspector.vertex());
    } else {
      Vertex left = inspector.vertex().left_child();
      Vertex right = inspector.vertex().right_child();

      rules_.emplace_back<std::initializer_list<RuleLetter>>({
        get_letter(left),
        get_letter(right)
      });
      rules_.back().debug_id = rules_.size() - 1;

      vertex_rules_.insert(
        std::make_pair(
          inspector.vertex(),
          &(rules_.back())
        )
      );
    }
    inspector.next();
  }
}

size_t JezRules::remove_crossing_blocks() {
  size_t eliminated_blocks_count = 0;
  for (auto& rule : rules_) {
    auto current = rule.begin();

    while (current != rule.end()) {
      assert(!current->is_empty_nonterminal());
      assert(current->is_valid());
      assert(!current->is_empty_terminal());
      auto next = std::next(current);
      assert(next == rule.end() || !next->is_empty_nonterminal());
      assert(next == rule.end() || !next->is_empty_terminal());

      if (next != rule.end() &&
          current->last_terminal_letter_id() ==
          next->first_terminal_letter_id()) {

        current = rule.pop_last_from_letter(current);

        assert(!current->is_nonterminal());

        next = std::next(current);

        assert(next != rule.end());
        assert(current != next);
        assert(current->last_terminal_letter_id() ==
               next->first_terminal_letter_id());

        next = rule.pop_first_from_letter(next);

        assert(current->is_valid());
        assert(next->is_valid());
        assert(!current->is_empty_terminal() && !current->is_nonterminal());
        assert(!next->is_empty_terminal() && !next->is_nonterminal());
        assert(current->last_terminal_letter_id() == next->first_terminal_letter_id());
        current = rule.merge_letters(current, next);

        if (current->is_empty_terminal()) {
          ++eliminated_blocks_count;
          current = rule.remove_empty_letter(current).second;
        }

        if (rule.empty()) {
          for (auto & occurence : rule.nonterminal_index_) {
            occurence.rule_->remove_empty_letter(occurence.letter_);
          }
        }
      } else {
        current = next;
      }
    }
  }
  return eliminated_blocks_count;
}

std::vector<LetterPosition> JezRules::list_blocks() {
  std::vector<LetterPosition> blocks;

  for (auto& rule : rules_) {
    for (auto letter = rule.begin(); letter != rule.end(); ++letter) {
      assert(!letter->is_empty_nonterminal() && !letter->is_empty_terminal());
      if (letter->is_power()) {
        blocks.push_back(LetterPosition(&rule, letter));
      }
    }
  }

  std::stable_sort(
    blocks.begin(),
    blocks.end(),
    [](const LetterPosition& first, const LetterPosition& second) -> bool {
      return first.letter_->terminal_id() < second.letter_->terminal_id() ||
          (first.letter_->terminal_id() == second.letter_->terminal_id() &&
           mad_sorts::reverse_bit_mpz_less(first.letter_->terminal_power(), second.letter_->terminal_power()));
    }
  );

  return blocks;
}

void JezRules::compress_blocks(const std::vector<LetterPosition>& blocks) {
  //LetterPower last_power = 0;
  TerminalId last_id = 0;

  Vertex current_terminal_vertex;
  std::vector<std::pair<LetterPower, Vertex*>> current_terminal_powers_;

  auto block = blocks.begin();
  while (block != blocks.end()) {
    assert(!block->letter_->is_nonterminal());
    if (!block->letter_->is_valid()) {
      ++block;
      continue;
    }

    assert(!block->rule_->empty());
    assert(block->letter_->is_power());

    if (last_id != block->letter_->terminal_id()) {
      last_id = block->letter_->terminal_id();
      current_terminal_vertex = terminal_vertices_[last_id];
      assert(current_terminal_powers_.empty());
    }

    if (current_terminal_powers_.empty() ||
        mpz_cmpabs(current_terminal_powers_.back().first.get_mpz_t(),
                   block->letter_->terminal_power().get_mpz_t()) != 0
    ) {
      this->terminal_vertices_.insert(
        this->terminal_vertices_.end(),
        std::make_pair(
          next_fresh_terminal(),
          Vertex()
        )
      );

      current_terminal_powers_.push_back(
        std::make_pair(
          abs(block->letter_->terminal_power()),
          &(std::prev(this->terminal_vertices_.end())->second)
        )
      );
    }

    block->rule_->compress_power(
        block->letter_,
        last_terminal()
    );

    ++block;

    if (block == blocks.end() ||
        last_id != block->letter_->terminal_id()
    ) {
      bool continue_iterations = true;
      while (continue_iterations) {
        Vertex last_vertex;
        Vertex last_power = current_terminal_vertex;
        continue_iterations = false;
        for (auto& terminal_power : current_terminal_powers_) {
          if ((terminal_power.first & 1) != 0) {
            if (*terminal_power.second != last_vertex) {
              last_vertex = *terminal_power.second;
              last_power = NonterminalVertex(last_vertex, current_terminal_vertex);
            }
            *terminal_power.second = last_power;
          }

          terminal_power.first >>= 1;

          if (terminal_power.first != 0) {
            continue_iterations = true;
          }
        }
        current_terminal_vertex = NonterminalVertex(
          current_terminal_vertex,
          current_terminal_vertex
        );
      }
      current_terminal_powers_.clear();
    }
  }
}

//void JezRules::empty_cleanup() {
//  for (auto& rule : rules_) {
//    for (auto current = rule.begin(); current != rule.end(); ) {
//      auto next = std::next(current);
//      if (current->is_empty_nonterminal()) {
//        rule.remove_empty_letter(current);
//      }
//      current = next;
//    }
//  }
//}

OneStepPairs::OneStepPairs(JezRules* rules)
  : rules_(rules)
{
  std::vector<std::tuple<
      TerminalId, //first letter
      TerminalId, //second letter
      LetterPosition>> all_pairs; //reference to position

  for (auto& rule : rules->rules_) {
    for (
        auto current = rule.begin(), next = std::next(rule.begin());
        next != rule.end();
        current = next, ++next
    ) {
      if (
          current->last_terminal_letter_id() != next->first_terminal_letter_id()
          //&& !current->is_power() && !next->is_power()
      ) {
        if (current->last_terminal_letter_sign() > 0 && next->first_terminal_letter_sign() > 0) { //12; 21
          all_pairs.emplace_back(
              current->last_terminal_letter_id(),
              next->first_terminal_letter_id(),
              LetterPosition(&rule, current)
          );
        } else if (current->last_terminal_letter_sign() < 0 && next->first_terminal_letter_sign() < 0) {//-1-2 -> 21, -2-1 -> 12
          all_pairs.emplace_back(
              next->first_terminal_letter_id(),
              current->last_terminal_letter_id(),
              LetterPosition(&rule, current)
          );
        } else if (current->last_terminal_letter_sign() < 0) {
          all_pairs.emplace_back(
              -std::max(current->last_terminal_letter_id(), next->first_terminal_letter_id()), //-21; -12->-21
              std::min(current->last_terminal_letter_id(), next->first_terminal_letter_id()),
              LetterPosition(&rule, current)
          );
        } else {
          all_pairs.emplace_back(
              std::min(current->last_terminal_letter_id(), next->first_terminal_letter_id()), //1-2, 2-1->1-2
              -std::max(current->last_terminal_letter_id(), next->first_terminal_letter_id()),
              LetterPosition(&rule, current)
          );
        }

      }
    }
  }

  std::stable_sort(
      all_pairs.begin(),
      all_pairs.end(),
      [] (const std::tuple<TerminalId, TerminalId, LetterPosition>& first,
          const std::tuple<TerminalId, TerminalId, LetterPosition>& second) {
      if (mad_sorts::signed_id_less(std::get<0>(first), std::get<0>(second))) {
        return true;
      } else if (std::get<0>(first) == std::get<0>(second)) {
        return mad_sorts::signed_id_less(std::get<1>(first), std::get<1>(second));
      }
      return false;
    }
  );

  if (all_pairs.empty()) {
    return;
  }

  auto left_letter_iterator = pairs_.end();
  auto right_letter_iterator = pairs_.end();

  decltype(right_letter_iterator->left_letters_.end()) left_list_current_letter;
  decltype(left_letter_iterator->right_letters_.end()) right_list_current_letter;

  for (const auto& pair: all_pairs) {
    while (
        left_letter_iterator != pairs_.end() &&
        mad_sorts::signed_id_less(left_letter_iterator->id_, std::get<0>(pair))
    ) {
      ++left_letter_iterator;
      right_letter_iterator = pairs_.begin();
      right_list_current_letter = left_letter_iterator->right_letters_.begin();
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    if (left_letter_iterator == pairs_.end() ||
        left_letter_iterator->id_ != std::get<0>(pair)
    ) {
      left_letter_iterator = pairs_.emplace(left_letter_iterator, std::get<0>(pair));
      right_letter_iterator = pairs_.begin();
      right_list_current_letter = left_letter_iterator->right_letters_.begin();
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    while (
        right_list_current_letter != left_letter_iterator->right_letters_.end() &&
        mad_sorts::signed_id_less(right_list_current_letter->id_, std::get<1>(pair))
    ) {
      ++right_list_current_letter;
    }

    if (
        right_list_current_letter == left_letter_iterator->right_letters_.end() ||
        right_list_current_letter->id_ != std::get<1>(pair)
    ) {
      right_list_current_letter = left_letter_iterator->right_letters_.insert(
          right_list_current_letter,
          std::get<1>(pair)
      );
    }

    right_list_current_letter->occurencies.push_back(std::get<2>(pair));

    while (
        right_letter_iterator != pairs_.end() &&
        mad_sorts::signed_id_less(right_letter_iterator->id_, std::get<1>(pair))
    ) {
      ++right_letter_iterator;
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    if (right_letter_iterator == pairs_.end() ||
        right_letter_iterator->id_ != std::get<1>(pair)
    ) {
      right_letter_iterator = pairs_.insert(
          right_letter_iterator,
          std::get<1>(pair)
      );

      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    while (
        left_list_current_letter != right_letter_iterator->left_letters_.end() &&
        mad_sorts::signed_id_less(*left_list_current_letter, std::get<0>(pair))
    ) {
      ++left_list_current_letter;
    }

    if (
        left_list_current_letter == right_letter_iterator->left_letters_.end() ||
        *left_list_current_letter != std::get<0>(pair)
    ) {
      left_list_current_letter = right_letter_iterator->left_letters_.insert(
          left_list_current_letter,
          std::get<0>(pair)
      );
    }
  }
}


OneStepPairs::GreedyLettersSeparation::GreedyLettersSeparation(const OneStepPairs& pairs)
  : flipped_(false)
  , id_shift_(0)
  , max_id_(0)
{
  const std::list<LetterLeftRight>& pairs_ = pairs.pairs_;
  auto letter_pairs = pairs_.begin();

  if (!pairs_.empty()) {
    id_shift_ = abs(pairs_.front().id_);
    max_id_ = abs(pairs_.back().id_);
  }

  size_t new_pairs_if_right = 0, new_pairs_if_left = 0;
  while (letter_pairs != pairs_.end()) {
    TerminalId letter = letter_pairs->id_;

    for (const auto& right_letter : letter_pairs->right_letters_) {
      assert(right_letter.id_ != letter);
      if (abs(right_letter.id_) > abs(letter) || is_right_letter(right_letter.id_)) {
        if (letter > 0) {
          ++new_pairs_if_left;
        } else {
          ++new_pairs_if_right;
        }
      }
    }

    for (const auto& left_letter : letter_pairs->left_letters_) {
      assert(left_letter != letter);
      if (abs(left_letter) > abs(letter) || is_left_letter(left_letter)) {
        if (letter > 0) {
          ++new_pairs_if_right;
        } else {
          ++new_pairs_if_left;
        }
      }
    }

    ++letter_pairs;
    if ((letter_pairs == pairs_.end() || abs(letter_pairs->id_) != abs(letter)) && (new_pairs_if_left > 0 || new_pairs_if_right > 0)) {
      right_letters_.resize(max_id_ - id_shift_ + 1, 0);
      if (new_pairs_if_left < new_pairs_if_right) {
        right_letters_[abs(letter) - id_shift_] = 1;
      }
      new_pairs_if_left = 0;
      new_pairs_if_right = 0;
    }
    if (empty() && letter_pairs != pairs_.end() && abs(letter_pairs->id_) != abs(letter)) {
      id_shift_ = abs(letter_pairs->id_);
    }
  }
}

TerminalId OneStepPairs::compress_pair(
    TerminalId first,
    TerminalId second,
    const std::vector<LetterPosition>& occurencies
) {
  TerminalId pair_terminal_id = rules_->next_fresh_terminal();

  for (auto& occurence : occurencies) {
    auto first_letter = occurence.letter_;
    auto second_letter = std::next(first_letter);

    assert(first_letter->is_valid());
    assert(!occurence.rule_->empty());
    assert(!first_letter->is_nonterminal());
    assert(!first_letter->is_power());
    assert(second_letter != occurence.rule_->end());
    assert(!second_letter->is_nonterminal());
    assert(!second_letter->is_power());
    assert((first_letter->terminal_sgn() * first > 0 && second_letter->terminal_sgn() * second > 0) || (first_letter->terminal_sgn() * first < 0 && second_letter->terminal_sgn() * second < 0));

    occurence.rule_->compress_pair(
        first_letter,
        second_letter,
        first_letter->terminal_id() == abs(first) ? pair_terminal_id : -pair_terminal_id
    );
  }
  return pair_terminal_id;
}

void OneStepPairs::remove_crossing(
    const OneStepPairs::GreedyLettersSeparation& letters_separation
) {
  std::vector<std::tuple<TerminalId, TerminalId, LetterPosition>> occurencies;

  for (auto & rule : rules_->rules_) {
    if (rule.empty()) {
      continue;
    }
    auto current = rule.begin();

    assert(!current->is_empty_nonterminal());

    auto next = std::next(current);
    while (current != rule.end() && next != rule.end()) {
      assert(!next->is_empty_nonterminal());

      if (letters_separation.is_left_letter(current->last_terminal_letter_id() * current->last_terminal_letter_sign()) &&
          letters_separation.is_right_letter(next->first_terminal_letter_id() * next->first_terminal_letter_sign()))
      {
        current = rule.pop_last_from_letter(current);
        assert(current->is_valid());
        assert(!current->is_empty_nonterminal());
        assert(!current->is_empty_terminal());

        next = std::next(current);

        assert(next != rule.end());
        next = rule.pop_first_from_letter(next);

        assert(next->is_valid());
        assert(!next->is_empty_nonterminal());
        assert(!next->is_empty_terminal());

        assert(letters_separation.is_left_letter(current->last_terminal_letter_id() * current->last_terminal_letter_sign()));
        assert(letters_separation.is_right_letter(next->first_terminal_letter_id() * next->first_terminal_letter_sign()));

        assert(std::next(current) == next);

        if (current->last_terminal_letter_sign() > 0 && next->first_terminal_letter_sign() > 0) { //12; 21
          occurencies.emplace_back(
              current->last_terminal_letter_id(),
              next->first_terminal_letter_id(),
              LetterPosition(&rule, current)
          );
        } else if (current->last_terminal_letter_sign() < 0 && next->first_terminal_letter_sign() < 0) {//-1-2 -> 21, -2-1 -> 12
          occurencies.emplace_back(
              next->first_terminal_letter_id(),
              current->last_terminal_letter_id(),
              LetterPosition(&rule, current)
          );
        } else if (current->last_terminal_letter_sign() < 0) {
          occurencies.emplace_back(
              -std::max(current->last_terminal_letter_id(), next->first_terminal_letter_id()), //-21; -12->-21
              std::min(current->last_terminal_letter_id(), next->first_terminal_letter_id()),
              LetterPosition(&rule, current)
          );
        } else {
          occurencies.emplace_back(
              std::min(current->last_terminal_letter_id(), next->first_terminal_letter_id()), //1-2, 2-1->1-2
              -std::max(current->last_terminal_letter_id(), next->first_terminal_letter_id()),
              LetterPosition(&rule, current)
          );
        }
      }
      current = next;
      next = std::next(current);
    }
  }

  std::stable_sort(
      occurencies.begin(),
      occurencies.end(),
      [] (const std::tuple<TerminalId, TerminalId, LetterPosition>& first,
          const std::tuple<TerminalId, TerminalId, LetterPosition>& second) {
        if (mad_sorts::signed_id_less(std::get<0>(first), std::get<0>(second))) {
          return true;
        } else if (std::get<0>(first) == std::get<0>(second)) {
          return mad_sorts::signed_id_less(std::get<1>(first), std::get<1>(second));
        }
        return false;
      }
  );

  auto current_pair = occurencies.begin();

  for (auto& pair : pairs_) {
    if (current_pair != occurencies.end() &&
        !mad_sorts::signed_id_less(pair.id_, std::get<0>(*current_pair))) {
      assert(pair.id_ == std::get<0>(*current_pair));
      for (auto & right_letter : pair.right_letters_) {
        if (letters_separation.is_right_letter(right_letter.id_)) {
          right_letter.occurencies.clear();
          if (current_pair != occurencies.end() &&
              !mad_sorts::signed_id_less(pair.id_, std::get<0>(*current_pair)) &&
              !mad_sorts::signed_id_less(right_letter.id_, std::get<1>(*current_pair))) {
            assert(right_letter.id_ == std::get<1>(*current_pair));
            while (pair.id_ == std::get<0>(*current_pair) &&
                right_letter.id_ == std::get<1>(*current_pair)) {
              right_letter.occurencies.push_back(std::get<2>(*current_pair));
              ++current_pair;

              assert(!right_letter.occurencies.back().rule_->empty());
              assert(!right_letter.occurencies.back().letter_->is_nonterminal());
              assert(!right_letter.occurencies.back().letter_->is_power());
              assert(std::next(right_letter.occurencies.back().letter_) !=
                  right_letter.occurencies.back().rule_->begin());
              assert(!std::next(right_letter.occurencies.back().letter_)->is_nonterminal());
              assert(!std::next(right_letter.occurencies.back().letter_)->is_power());

              if (current_pair == occurencies.end()) {
                break;
              }
            }
          }
        }
      }
    } else {
      if (letters_separation.is_left_letter(pair.id_)) {
        for (auto & right_letter : pair.right_letters_) {
          if (letters_separation.is_right_letter(right_letter.id_)) {
            right_letter.occurencies.clear();
          }
        }
      }
    }
  }
}


void OneStepPairs::compress_pairs_from_letter_lists(
    const OneStepPairs::GreedyLettersSeparation& letters_separation
) {
  for (auto& pair : pairs_) {
    if (letters_separation.is_left_letter(pair.id_)) {
      auto right_letter = pair.right_letters_.begin();

      while (right_letter != pair.right_letters_.end()) {
        if (letters_separation.is_right_letter(right_letter->id_)) {
          auto terminal = compress_pair(
              pair.id_,
              right_letter->id_,
              right_letter->occurencies
          );

          rules_->terminal_vertices_.insert(std::make_pair(
              terminal,
              NonterminalVertex(
                  pair.id_ > 0 ? rules_->terminal_vertices_[pair.id_] : rules_->terminal_vertices_[-pair.id_].negate(),
                  right_letter->id_ > 0 ? rules_->terminal_vertices_[right_letter->id_] : rules_->terminal_vertices_[-right_letter->id_].negate()
              )
          ));

          right_letter = pair.right_letters_.erase(right_letter);
        } else {
          ++right_letter;
        }
      }
    }

    if (letters_separation.is_right_letter(pair.id_)) {
      auto left_letter = pair.left_letters_.begin();

      while (left_letter != pair.left_letters_.end()) {
        if (letters_separation.is_left_letter(*left_letter)) {
          left_letter = pair.left_letters_.erase(left_letter);
        } else {
          ++left_letter;
        }
      }
    }
  }
}

Vertex JezRules::get_exponent(Vertex vertex, LetterPower power) {
  Vertex result;

  assert(power > 0);

  while (power != 0) {
    if ((power & 1) != 0) {
      if (result) {
        result = NonterminalVertex(result, vertex);
      } else {
        result = vertex;
      }
    }
    vertex = NonterminalVertex(vertex, vertex);
    power >>= 1;
  }

  return result;
}

void RuleLetter::debug_print(std::ostream* os) const {
  if (this->is_nonterminal()) {
    (*os) << "(" << this->nonterminal_rule()->debug_id << ")";
  } else {
    (*os) << this->terminal_id() << "^" << this->terminal_power();
  }
}

void Rule::debug_print(::std::ostream* os) const {
  for (auto& letter : letters_) {
    letter.debug_print(os);
  }
}

void Rule::debug_print_exposed(::std::ostream* os) const {
  for (auto& letter : letters_) {
    if (letter.is_nonterminal()) {
      letter.nonterminal_rule_->debug_print_exposed(os);
    } else {
      for (LetterPower i = 0; i < letter.terminal_power(); ++i) {
        (*os) << letter.terminal_id() << ',';
      }
    }
  }
}

void JezRules::debug_print(std::ostream* os) const {
  for (auto& rule : rules_) {
    (*os) << rule.debug_id << ": ";
    (*os) << rule.first_terminal_id() << ".." <<
        rule.last_terminal_id() << ": ";
    for (auto& letter : rule) {
      letter.debug_print(os);
      (*os) << ", ";
    }
    (*os) << std::endl;
  }
}

void OneStepPairs::GreedyLettersSeparation::debug_print(std::ostream* os) const {
  (*os) << "\nGreedyPairs:\nshift: " << id_shift_ << "\nLeft: ";

  for(size_t i = 0; i < right_letters_.size(); ++i) {
    if (!right_letters_[i]) {
      (*os) << i + id_shift_<< ',';
    }
  }

  (*os) << "\nRight: ";
  for(size_t i = 0; i < right_letters_.size(); ++i) {
    if (right_letters_[i]) {
      (*os) << i + id_shift_ << ',';
    }
  }
  (*os) << std::endl;
}

Vertex normal_form(Vertex root) {
  if (root.height() < 2) {
    return root;
  }
  auto rules = JezRules::create(root);
  Rule& root_rule = *(rules->vertex_rules_[root]);

  while (root_rule.size() > 1 || (
           root_rule.size() == 1 && (
             root_rule.begin()->is_power() ||
             root_rule.begin()->is_nonterminal()
           )
        )) {

//    std::cout << "\n=================\n\nCurrent rules:" << std::endl;
//    rules.debug_print(&std::cout);

    rules->remove_crossing_blocks();
    OneStepPairs pairs(rules.get());
//
//    std::cout << "Rules after RemCrBlocks: " << std::endl;
//    rules.debug_print(&std::cout);

    auto blocks = rules->list_blocks();
//    std::cout << "\nFound blocks: " << std::endl;
//    for (auto& block : blocks) {
//      std::cout << block.rule_->debug_id << ':';
//      block.letter_->debug_print(&std::cout);
//      std::cout << std::endl;
//    }

    rules->compress_blocks(blocks);

//    std::cout << "Rules after CompressBlocks: " << std::endl;
//    rules.debug_print(&std::cout);

    OneStepPairs::GreedyLettersSeparation letters_separation(pairs);
//    letters_separation.debug_print(&std::cout);

    while (!letters_separation.empty()) {
      pairs.remove_crossing(letters_separation);
      pairs.compress_pairs_from_letter_lists(letters_separation);
//      std::cout << "Rules after first compression: " << std::endl;
//      rules.debug_print(&std::cout);

      letters_separation.flip();
      pairs.remove_crossing(letters_separation);
      pairs.compress_pairs_from_letter_lists(letters_separation);
//      std::cout << "Rules after second compression: " << std::endl;
//      rules.debug_print(&std::cout);
      letters_separation = OneStepPairs::GreedyLettersSeparation(pairs);
//      letters_separation.debug_print(&std::cout);
    }

    //rules->empty_cleanup();
  }

  Rule::collect_garbage();
  return rules->terminal_vertices_[root_rule.first_terminal_id()];
}

//bool JezReducingRules::remove_empty_terminals() {
//  bool reduction_possible = false;
//  for (auto& rule : rules_) {
//    for (auto current = rule.begin(); current != rule.end(); ) {
//      auto next = std::next(current);
//      assert(!current->is_empty_nonterminal());
//      if (next != rule.end() &&
//          current->last_terminal_letter_id() == next->first_terminal_letter_id() &&
//          current->last_terminal_letter_sign() != next->first_terminal_letter_sign()) {
//        reduction_possible = true;
//      }
//      if (current->is_empty_terminal()) {
//        current = rule.remove_empty_terminal(current).first;
//        reduction_possible = true;
//      } else {
//        current = next;
//      }
//    }
//  }
//  return reduction_possible;
//}

Vertex reduce(Vertex root) {
  if (root.height() < 2) {
    return root;
  }

  bool reduced = false;
  bool normalized = false;

  while (root && !reduced) {
    auto rules = JezReducingRules::create(root);

    Rule* root_rule = rules->vertex_rules_[root];

    while (!root_rule->is_trivial()) {
      std::cout << "\n=================\n\nCurrent rules:" << std::endl;
      rules->debug_print(&std::cout);
      reduced = !rules->remove_crossing_blocks();

      std::cout << "Rules after RemCrBlocks: " << std::endl;
      rules->debug_print(&std::cout);

      OneStepPairs pairs(rules.get());

      auto blocks = rules->list_blocks();
      std::cout << "\nFound blocks: " << std::endl;
      for (auto& block : blocks) {
        std::cout << block.rule_->debug_id << ':';
        block.letter_->debug_print(&std::cout);
        std::cout << std::endl;
      }

      rules->compress_blocks(blocks);

      std::cout << "Rules after CompressBlocks: " << std::endl;
      rules->debug_print(&std::cout);

      OneStepPairs::GreedyLettersSeparation letters_separation(pairs);
      letters_separation.debug_print(&std::cout);

      while (!letters_separation.empty()) {
        pairs.remove_crossing(letters_separation);
        pairs.compress_pairs_from_letter_lists(letters_separation);
        std::cout << "Rules after first compression: " << std::endl;
        rules->debug_print(&std::cout);

        letters_separation.flip();
        pairs.remove_crossing(letters_separation);
        pairs.compress_pairs_from_letter_lists(letters_separation);
        std::cout << "Rules after second compression: " << std::endl;
        rules->debug_print(&std::cout);
        letters_separation = OneStepPairs::GreedyLettersSeparation(pairs);
        letters_separation.debug_print(&std::cout);
      }
//      rules->empty_cleanup();
    }

    normalized = true;

    Rule::collect_garbage();
    if (!root_rule->empty()) {
      while (root_rule->begin()->is_nonterminal()) {
        root_rule = root_rule->begin()->nonterminal_rule();
      }

      root = rules->terminal_vertices_[root_rule->first_terminal_id()];

      if (root_rule->begin()->first_terminal_letter_sign() < 0) {
        root = root.negate();
      }
    } else {
      root = Vertex();
    }
  }
  return root;
}

}
}
}



