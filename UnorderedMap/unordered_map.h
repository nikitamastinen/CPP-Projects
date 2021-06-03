#include <vector>
#include <type_traits>
#include <stdexcept>
#include <iostream>

template<typename T, typename Allocator = std::allocator<T>>
class List {
  class Node {
  public:
    Node* next;
    Node* prev;
    T value;

    explicit Node() : next(nullptr), prev(nullptr), value(T()) {}
    explicit Node(const T& value) : next(nullptr), prev(nullptr), value(value) {}
  };

  Node* head;
  size_t length;
  using NAllocator = typename Allocator::template rebind<Node>::other;
  NAllocator allocator;
  Allocator t_allocator;
public:
  Node* insert(Node* node, const T& value) {
    Node* ins = std::allocator_traits<NAllocator>::allocate(allocator, 1);
    std::allocator_traits<NAllocator>::construct(allocator, ins, value);
    ins->next = node->next;
    ins->prev = node;
    node->next->prev = ins;
    node->next = ins;
    ++length;
    return ins;
  }

  Node* erase(Node* node) {
    if (node == head) {
      return head;
    }
    node->prev->next = node->next;
    node->next->prev = node->prev;
    Node* result = node->next;
    std::allocator_traits<NAllocator>::destroy(allocator, node);
    std::allocator_traits<NAllocator>::deallocate(allocator, node, 1);
    --length;
    return result;
  }

  void no_allocator_swap(List& other) {
    std::swap(length, other.length);
    std::swap(head, other.head);
  }
public:
  explicit List(const Allocator& t_allocator = Allocator()):
      length(0), t_allocator(t_allocator) {
    typename Allocator::template rebind<Node>::other n_allocator;
    allocator = n_allocator;
    head = std::allocator_traits<NAllocator>::allocate(allocator, 1);
    head->next = head;
    head->prev = head;
  }

  explicit List(
      size_t count,
      const T& value,
      const Allocator& t_allocator = Allocator()):
      length(count), t_allocator(t_allocator) {
    typename Allocator::template rebind<Node>::other n_allocator;
    allocator = n_allocator;
    head = std::allocator_traits<NAllocator>::allocate(allocator, 1);
    head->next = head;
    head->prev = head;
    Node* copy_before = head;
    Node* copy_after = head;
    for (size_t i = 0; i < count; i++) {
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator, 1);
      std::allocator_traits<NAllocator>::construct(allocator, copy_after, value);
      copy_before->next = copy_after;
      copy_after->prev = copy_before;
      copy_before = copy_after;
    }
    copy_after->next = head;
    head->prev = copy_after;
  }

  explicit List(
      size_t count,
      const Allocator& t_allocator = Allocator()):
      length(count), t_allocator(t_allocator) {
    typename Allocator::template rebind<Node>::other n_allocator;
    allocator = n_allocator;
    head = std::allocator_traits<NAllocator>::allocate(allocator, 1);
    head->next = head;
    head->prev = head;
    Node* copy_before = head;
    Node* copy_after = head;
    for (size_t i = 0; i < count; i++) {
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator, 1);
      std::allocator_traits<NAllocator>::construct(allocator, copy_after);
      copy_before->next = copy_after;
      copy_after->prev = copy_before;
      copy_before = copy_after;
    }
    copy_after->next = head;
    head->prev = copy_after;
  }

  List(const List& other) : length(other.length) {
    t_allocator = std::allocator_traits<Allocator>::
    select_on_container_copy_construction(other.t_allocator);
    allocator = std::allocator_traits<NAllocator>::
    select_on_container_copy_construction(other.allocator);
    head = std::allocator_traits<NAllocator>::allocate(allocator, 1);
    head->next = head;
    head->prev = head;
    Node* copy_before = head;
    Node* copy_after = head;
    Node* other_head = other.head;
    for (size_t i = 0; i < length; i++) {
      other_head = other_head->next;
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator, 1);
      std::allocator_traits<NAllocator>::
      construct(allocator, copy_after, other_head->value);
      copy_before->next = copy_after;
      copy_after->prev = copy_before;
      copy_before = copy_after;
    }
    copy_after->next = head;
    head->prev = copy_after;
  }

  void clear() {
    while (begin() != end()) {
      erase(begin());
    }
  }

  List(List&& other) noexcept : length(other.length) {
      t_allocator = std::move(other.t_allocator);
      allocator = std::move(other.allocator);
      auto new_head = std::allocator_traits<NAllocator>::allocate(allocator, 1);
      new_head->next = new_head;
      new_head->prev = new_head;
      other.length = 0;
      head = other.head;
      other.head = new_head;
  }

  List& operator=(const List& other) {
    if (this == &other) {
      return *this;
    }
    List copied = other;
    if (std::allocator_traits<Allocator>::
    propagate_on_container_copy_assignment::value) {
      t_allocator = other.t_allocator;
      allocator = other.allocator;
    }
    no_allocator_swap(copied);
    return *this;
  }

  List& operator=(List&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    List copied = std::move(other);
    if (std::allocator_traits<Allocator>::
    propagate_on_container_move_assignment::value) {
      t_allocator = std::move(other.t_allocator);
      allocator = std::move(other.allocator);
    }
    no_allocator_swap(copied);
    return *this;
  }

  ~List() {
    Node* copy_before = head->next;
    for (size_t i = 0; i < length; i++) {
      auto copy_after = copy_before->next;
      std::allocator_traits<NAllocator>::destroy(allocator, copy_before);
      std::allocator_traits<NAllocator>::deallocate(allocator, copy_before, 1);
      copy_before = copy_after;
    }
    std::allocator_traits<NAllocator>::deallocate(allocator, head, 1);
  }

  Allocator get_allocator() const {
    return t_allocator;
  }

  template<bool IsConst>
  class iterator_impl {
  public:
    std::conditional_t<IsConst, const Node*, Node*> it;
    friend class List;
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = T;
    using pointer = std::conditional_t<IsConst, const T*, T*>;;
    using reference = std::conditional_t<IsConst, const T&, T&>;
    using const_reference = const T&;
    using difference_type = std::ptrdiff_t;

    iterator_impl(Node* node) {
      it = node;
    }

    operator iterator_impl<true>() {
      return iterator_impl<true>(it);
    }

    std::conditional_t<IsConst, const_reference, reference>  operator*() const {
      return it->value;
    }

    std::conditional_t<IsConst, const pointer , pointer> operator->() const {
      return &(it->value);
    }

    iterator_impl& operator++() {
      it = it->next;
      return *this;
    }

    iterator_impl operator++(int) {
      auto copied = *this;
      it = it->next;
      return copied;
    }

    iterator_impl& operator--() {
      it = it->prev;
      return *this;
    }

    iterator_impl operator--(int) {
      auto copied = *this;
      it = it->prev();
      return copied;
    }

    bool operator==(const iterator_impl& other) const {
      return it == other.it;
    }

    bool operator!=(const iterator_impl& other) const {
      return !(*this == other);
    }
  };

  using const_iterator = iterator_impl<true>;
  using iterator = iterator_impl<false>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  size_t size() const {
    return length;
  }

  iterator begin() const {
    iterator it(head->next);
    return it;
  }

  const_iterator cbegin() const {
    const_iterator it(head->next);
    return it;
  }

  reverse_iterator rbegin() const {
    reverse_iterator it(head);
    return it;
  }

  const_reverse_iterator crbegin() const {
    const_reverse_iterator it(head);
    return it;
  }

  iterator end() const {
    iterator it(head);
    return it;
  }

  const_iterator cend() const {
    const_iterator it(head->next);
    return it;
  }

  reverse_iterator rend() const {
    reverse_iterator it(head->next);
    return it;
  }

  const_reverse_iterator crend() const {
    const_reverse_iterator it(head->next);
    return it;
  }

  void push_front(const T& value) {
    insert(head, value);
  }

  void push_back(const T& value) {
    insert(head->prev, value);
  }

  void pop_front() {
    erase(head->next);
  }

  void pop_back() {
    erase(head->prev);
  }

  iterator insert(const_iterator it, const T& value) {
    return iterator(insert(const_cast<Node*>((--it).it), value));
  }

  iterator erase(const_iterator it) {
    return iterator(erase(const_cast<Node*>(it.it)));
  }
};

template<
    typename Key,
    typename Value,
    typename Hash = std::hash<Key>,
    typename Equal = std::equal_to<Key>,
    typename Allocator = std::allocator<std::pair<const Key, Value>>
>
class UnorderedMap {
public:
  using NodeType = std::pair<const Key, Value>;
private:
  static const size_t INITIAL_SIZE_OF_CONTAINER = 100;
  struct NodeTypeWithCashedHash {
    NodeType* ptr;
    size_t hash = 0;

    NodeType* get() {
      return ptr;
    }
  };
  using ListIterator = typename List<NodeTypeWithCashedHash*>::iterator;
  using ConstListIterator = typename List<NodeTypeWithCashedHash*>::const_iterator;

  template<bool IsConst>
  class iterator_impl {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = NodeType;
    using difference_type = size_t;
    using pointer = typename std::conditional_t<IsConst, const NodeType*, NodeType*>;
    using reference = typename std::conditional_t<IsConst, const NodeType&, NodeType&>;
    using list_iterator =
    typename std::conditional_t<IsConst, ConstListIterator, ListIterator>;

    list_iterator it;
    reference operator*() const {
      return *((*it)->ptr);
    }

    pointer operator->() const {
      return (*it)->ptr;
    }

    list_iterator get() {
      return it;
    }

    iterator_impl(const list_iterator& node) : it(node) {}

    operator iterator_impl<true>() {
      return iterator_impl<true>(it);
    }

    iterator_impl& operator++() {
      ++it;
      return *this;
    }

    iterator_impl& operator--() {
      --it;
      return *this;
    }

    iterator_impl operator++(int) {
      auto copied = *this;
      ++it;
      return copied;
    }

    bool operator==(const iterator_impl& other) const {
      return it == other.it;
    }

    bool operator!=(const iterator_impl& other) const {
      return it != other.it;
    }
  };

public:
  using iterator = iterator_impl<false>;
  using const_iterator = iterator_impl<true>;

private:
  List<NodeTypeWithCashedHash*> elements_;
  std::vector<ListIterator> hash_array_;
  Hash hash_function_;
  using NAllocator =
  typename std::allocator_traits<Allocator>::template rebind_alloc<NodeTypeWithCashedHash>;
  NAllocator t_alloc_;
  Equal equal_key_;

public:
  float current_max_load_factor = 0.75;

  UnorderedMap() : hash_array_(INITIAL_SIZE_OF_CONTAINER, elements_.end()) {}

  UnorderedMap(const UnorderedMap& other) :
      hash_array_({INITIAL_SIZE_OF_CONTAINER + other.size(), elements_.end()}),
      hash_function_(other.hash_function_),
      t_alloc_(
          std::allocator_traits<Allocator>::select_on_container_copy_construction(other.t_alloc_)
      ),
      equal_key_(other.equal_key_),
      current_max_load_factor(other.current_max_load_factor) {

    for (auto it : other) {
      insert(it);
    }
  }

  UnorderedMap(UnorderedMap&& other) :
      elements_(std::move(other.elements_)),
      hash_array_(std::move(other.hash_array_)),
      hash_function_(std::move(other.hash_function_)),
      t_alloc_(std::move(other.t_alloc_)),
      equal_key_(std::move(other.equal_key_)),
      current_max_load_factor(std::move(other.current_max_load_factor))
  {}

private:
  void clear_list_elements() {
    if (elements_.size() == 0) return;

    for (iterator it = elements_.begin(); it != elements_.end(); ++it) {
      std::allocator_traits<NAllocator>::destroy(t_alloc_, &(*it));
      Allocator alloc = typename std::allocator_traits<NAllocator>::template rebind_alloc<NodeType>();
      std::allocator_traits<Allocator>::deallocate(alloc, (*(it.get()))->ptr, 1);
      std::allocator_traits<NAllocator>::deallocate(t_alloc_, *it.get(), 1);
    }
    elements_.clear();
  }

  void swap_and_kill(UnorderedMap&& other) {
    elements_ = std::move(other.elements_);
    hash_array_ = std::move(other.hash_array_);
    hash_function_ = std::move(other.hash_function_);
    clear_list_elements();
    equal_key_ = std::move(other.equal_key_);
    current_max_load_factor = std::move(other.current_max_load_factor);
  }

public:
  UnorderedMap& operator=(const UnorderedMap& other) {
    if (this == &other) {
      return *this;
    }

    UnorderedMap copied = other;
    if (std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) {
      t_alloc_ = other.t_alloc_;
    }
    swap_and_kill(std::move(copied));
    return *this;
  }

  UnorderedMap& operator=(UnorderedMap&& other) noexcept {
    if (this == &other) {
      return *this;
    }

    if (std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value) {
      t_alloc_ = std::move(other.t_alloc_);
    }
    swap_and_kill(std::move(other));
    return *this;
  }

  ~UnorderedMap() {
    clear_list_elements();
  }

  size_t size() const {
    return elements_.size();
  }

  size_t bucket_id(const Key& key) const {
    return hash_function_(key) % hash_array_.size();
  }

  float load_factor() const {
    return static_cast<float>(elements_.size()) / hash_array_.size();
  }

  float load_factor_after_insert() const {
    return static_cast<float>(elements_.size() + 1) / hash_array_.size();
  }

  void max_load_factor(float value) {
    current_max_load_factor = value;
  }

  float max_load_factor() {
    return current_max_load_factor;
  }

private:
  void update() {
    if (load_factor_after_insert() > current_max_load_factor) {
      reserve(hash_array_.size() * 2 + 1);
    }
  }

public:
  void reserve(size_t count) {
    if (count > hash_array_.size()) {
      rehash(static_cast<size_t>(static_cast<float>(count) / max_load_factor() + 1));
    }
  }

  void rehash(size_t count) {
    hash_array_.clear();
    List<NodeTypeWithCashedHash*> copied = std::move(elements_);
    hash_array_.resize(count, elements_.end());

    for (ListIterator it = copied.begin(); it != copied.end(); ++it) {
      size_t hsh = bucket_id((*it)->ptr->first);
      ListIterator& elem = hash_array_[hsh];
      if (elem == elements_.end()) {
        elem = elements_.insert(elements_.cend(), *it);
      } else {
        (*elem)->hash = hsh;
        elem = elements_.insert(elem, *it);
      }
    }
  }

  template<class... Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    NodeTypeWithCashedHash* mover = std::allocator_traits<NAllocator>::allocate(t_alloc_, 1);
    Allocator alloc;
    auto new_node = std::allocator_traits<Allocator>::allocate(alloc, 1);
    std::allocator_traits<Allocator>::construct(alloc, new_node, std::forward<Args>(args)...);

    mover->ptr = new_node;
    mover->hash = bucket_id(mover->ptr->first);
    iterator result = find(mover->ptr->first);

    if (result != elements_.end()) {
      std::allocator_traits<NAllocator>::destroy(t_alloc_, mover);
      std::allocator_traits<Allocator>::deallocate(alloc, mover->ptr, 1);
      std::allocator_traits<NAllocator>::deallocate(t_alloc_, mover, 1);
      return {result, false};
    }

    update();

    ListIterator& elem = hash_array_[bucket_id(mover->ptr->first)];
    if (elem == elements_.end()) {
      elements_.insert(elements_.cend(), mover);
    } else {
      elements_.insert(elem, mover);
    }
    return {elem, true};
  }

  std::pair<iterator, bool> insert(const NodeType& value) {
    iterator result = find(value.first);
    if (result != elements_.end()) {
      return {result, false};
    }
    update();
    auto hsh = bucket_id(value.first);
    ListIterator& elem = hash_array_[hsh];
    NodeTypeWithCashedHash* copied = std::allocator_traits<NAllocator>::allocate(t_alloc_, 1);
    Allocator alloc;
    auto new_node = std::allocator_traits<Allocator>::allocate(alloc, 1);
    std::allocator_traits<Allocator>::construct(alloc, new_node, value);
    copied->ptr = new_node;
    copied->hash = bucket_id(value.first);
    if (elem == elements_.end()) {
      elem = elements_.insert(elements_.end(), copied);
    } else {
      elem = elements_.insert(elem, copied);
    }
    return {elem, true};
  }

  template<typename NodePair>
  std::pair<iterator, bool> insert(NodePair&& value) {
    iterator result = find(value.first);
    if (result != iterator(elements_.end())) {
      return {result, false};
    }

    update();

    ListIterator& elem = hash_array_[bucket_id(value.first)];
    NodeTypeWithCashedHash* mover = std::allocator_traits<NAllocator>::allocate(t_alloc_, 1);
    Allocator alloc;
    auto new_node = std::allocator_traits<Allocator>::allocate(alloc, 1);
    std::allocator_traits<Allocator>::construct(alloc, new_node, std::forward<NodePair>(value));
    mover->ptr = new_node;
    mover->hash = bucket_id(value.first);

    if (elem == elements_.end()) {
      elem = elements_.insert(elements_.end(), mover);
    } else {
      elem = elements_.insert(elem, mover);
    }
    return {elem, true};
  }

  template<typename Input>
  void insert(Input first, Input last) {
    for (; first != last; insert(*first++));
  }

  bool erase(const Key& key) {
    iterator it = find(key);
    if (it == elements_.end()) {
      return false;
    }

    erase(it);
    return true;
  }

  iterator erase(const_iterator it) {
    size_t ind = bucket_id(it->first);
    if (const_iterator(hash_array_[ind]) != it) {
      return elements_.erase(it.it);
    }

    auto nit = elements_.erase(it.it);
    hash_array_[ind] = (
        nit != elements_.end() && bucket_id((*nit)->ptr->first) == ind ? nit : elements_.end()
    );
    return nit;
  }

  iterator erase(const_iterator first, const_iterator last) {
    const_iterator stop = last;
    --stop;
    for (const_iterator it = first; it != last; it = erase(it)) {
      if (it == stop) {
        return erase(it);
      }
    }
    return end();
  }

  iterator find(const Key& key) {
    size_t hash = bucket_id(key);
    iterator it = hash_array_[hash];
    while (it != elements_.end() && (*it.get())->hash == hash) {
      if (equal_key_((*it.get())->ptr->first, key)) {
        return it;
      }
      ++it;
    }
    return elements_.end();
  }

  Value& at(const Key& key) {
    size_t hash = bucket_id(key);
    ListIterator it = hash_array_[hash];
    while (it != elements_.end() && (*it)->hash == hash) {
      if (equal_key_((*it)->ptr->first, key)) {
        return (*it)->ptr->second;
      }
      ++it;
    }
    throw std::out_of_range("Target element doesn't exists");
  }

  Value& operator[](const Key& key) {
    size_t hash = bucket_id(key);
    ListIterator it = hash_array_[hash];
    while (it != elements_.end() && (*it)->hash == hash) {
      if (equal_key_((*it)->ptr->first, key)) {
        return (*it)->ptr->second;
      }
      ++it;
    }

    iterator new_iter = insert({key, Value()}).first;
    return new_iter->second;
  }

  iterator begin() {
    return elements_.begin();
  }

  const_iterator begin() const {
    return elements_.cbegin();
  }

  iterator end() {
    return elements_.end();
  }

  const_iterator end() const {
    return elements_.cend();
  }

  const_iterator cbegin() const {
    return elements_.cbegin();
  }

  const_iterator cend() const {
    return elements_.cend();
  }
};