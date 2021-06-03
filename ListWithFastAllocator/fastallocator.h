//
// Created by Nikita Mastinen on 14.03.2021.
//

#pragma once

#include <type_traits>
#include <vector>
#include <memory>

template<size_t BlockSize>
class FixedAllocator {
  std::vector<char*> pool_items_;
  std::vector<char*> items_to_destroy_;
  static const size_t NUMBER_OF_OBJECTS_IN_BLOCK = 96;
public:
  static FixedAllocator<BlockSize>& instance() {
    static auto singleton = FixedAllocator<BlockSize>();
    return singleton;
  };

  void* allocate() {
    if (pool_items_.empty()) {
      char* new_pool = static_cast<char*>(
          ::operator new(BlockSize * NUMBER_OF_OBJECTS_IN_BLOCK)
      );
      items_to_destroy_.push_back(new_pool);
      for (size_t i = 0; i < NUMBER_OF_OBJECTS_IN_BLOCK; ++i) {
        pool_items_.push_back(new_pool + i * BlockSize);
      }
    }
    char* allocated = pool_items_.back();
    pool_items_.pop_back();
    return static_cast<void*>(allocated);
  }

  void deallocate(void* ptr) {
    pool_items_.push_back(static_cast<char*>(ptr));
  }

  ~FixedAllocator() {
    for (auto* pool_item : items_to_destroy_) {
      ::operator delete(pool_item);
    }
  }
};

template<typename T>
class FastAllocator {
public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  FastAllocator() = default;

  template<typename U>
  FastAllocator(const FastAllocator<U>&) {}

  ~FastAllocator() = default;

  T* allocate(size_t count) {
    if (count * sizeof(T) <= 8) {
      return static_cast<T*>(FixedAllocator<8>::instance().allocate());
    }
    if (count * sizeof(T) <= 24) {
      return static_cast<T*>(FixedAllocator<24>::instance().allocate());
    }
    if (count * sizeof(T) <= 32) {
      return static_cast<T*>(FixedAllocator<32>::instance().allocate());
    }
    return static_cast<T*>(::operator new(count * sizeof(T)));
  }

  void deallocate(T* ptr, size_t count) {
    if ((count * sizeof(T)) <= 8) {
      FixedAllocator<8>::instance().deallocate(ptr);
      return;
    }
    if ((count * sizeof(T)) <= 24) {
      FixedAllocator<24>::instance().deallocate(ptr);
      return;
    }
    if ((count * sizeof(T)) <= 32) {
      FixedAllocator<32>::instance().deallocate(ptr);
      return;
    }
    ::operator delete(ptr);
  }
};

template<typename U, typename V> bool operator!=(
    const FastAllocator<U>& x, const FastAllocator<V>& y
) {
  return !(x == y);
}

template<typename U, typename V> bool operator==(
    const FastAllocator<U>&, const FastAllocator<V>&
) {
  return true;
}

template<typename T, typename Allocator = std::allocator<T>>
class List {
private:
  class Node {
  public:
    Node* next;
    Node* prev;
    T value;

    explicit Node(): next(nullptr), prev(nullptr), value(T()) {}
    explicit Node(const T& value): next(nullptr), prev(nullptr), value(value) {}
  };

  Node* head_;
  size_t length_;

  using NAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
  NAllocator allocator_;

  Allocator object_allocator_;


  void link_(Node* before, Node* after) {
    before->next = after;
    after->prev = before;
  }

  void cycle_(Node* node) {
    node->next = node;
    node->prev = node;
  }

  void insert_(Node* node, const T& value) {
    Node* ins = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
    std::allocator_traits<NAllocator>::construct(allocator_, ins, value);
    ins->prev = node;
    link_(ins, node->next);
    node->next = ins;
    ++length_;
  }

  void erase_(Node* node) {
    link_(node->prev, node->next);
    std::allocator_traits<NAllocator>::destroy(allocator_, node);
    std::allocator_traits<NAllocator>::deallocate(allocator_, node, 1);
    --length_;
  }

  void swap_(List& other) {
    std::swap(length_, other.length_);
    std::swap(head_, other.head_);
  }
public:
  explicit List(const Allocator& object_allocator = Allocator()):
      length_(0), object_allocator_(object_allocator) {
    head_ = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
    cycle_(head_);
  }

  explicit List(
      size_t count,
      const T& value,
      const Allocator& object_allocator = Allocator()) : length_(count), object_allocator_(object_allocator) {
    head_ = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
    cycle_(head_);

    Node* copy_before = head_;
    Node* copy_after = head_;

    for (size_t i = 0; i < count; i++) {
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
      std::allocator_traits<NAllocator>::construct(allocator_, copy_after, value);
      link_(copy_before, copy_after);
      copy_before = copy_after;
    }

    link_(copy_after, head_);
  }

  explicit List(
      size_t count,
      const Allocator& object_allocator = Allocator()) : length_(count), object_allocator_(object_allocator) {
    head_ = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
    cycle_(head_);

    Node* copy_before = head_;
    Node* copy_after = head_;

    for (size_t i = 0; i < count; i++) {
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
      std::allocator_traits<NAllocator>::construct(allocator_, copy_after);
      link_(copy_before, copy_after);
      copy_before = copy_after;
    }

    link_(copy_after, head_);
  }

  List(const List& other): length_(other.length_) {
    object_allocator_ = std::allocator_traits<Allocator>::
    select_on_container_copy_construction(other.object_allocator_);
    allocator_ = std::allocator_traits<NAllocator>::
    select_on_container_copy_construction(other.allocator_);

    head_ = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
    cycle_(head_);

    Node* copy_before = head_;
    Node* copy_after = head_;
    Node* other_head = other.head_;

    for (size_t i = 0; i < length_; i++) {
      other_head = other_head->next;
      copy_after = std::allocator_traits<NAllocator>::allocate(allocator_, 1);
      std::allocator_traits<NAllocator>::
      construct(allocator_, copy_after, other_head->value);
      link_(copy_before, copy_after);
      copy_before = copy_after;
    }

    link_(copy_after, head_);
  }

  List& operator=(const List& other) {
    if (this == &other) {
      return *this;
    }

    List copied = other;
    if (std::allocator_traits<Allocator>::
    propagate_on_container_copy_assignment::value) {
      object_allocator_ = other.object_allocator_;
      allocator_ = other.allocator_;
    }

    swap_(copied);
    return *this;
  }

  ~List() {
    Node* copy_before = head_->next;
    for (size_t i = 0; i < length_; i++) {
      auto* copy_after = copy_before->next;
      std::allocator_traits<NAllocator>::destroy(allocator_, copy_before);
      std::allocator_traits<NAllocator>::deallocate(allocator_, copy_before, 1);
      copy_before = copy_after;
    }
    std::allocator_traits<NAllocator>::deallocate(allocator_, head_, 1);
  }

  Allocator get_allocator() const {
    return object_allocator_;
  }

  template<bool IsConst>
  class iterator_impl {
  private:
    std::conditional_t<IsConst, const Node*, Node*> it;
    friend class List;

  public:
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

    reference operator*() const {
      return it->value;
    }

    pointer operator->() const {
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
    return length_;
  }

  iterator begin() const {
    iterator it(head_->next);
    return it;
  }

  const_iterator cbegin() const {
    const_iterator it(head_->next);
    return it;
  }

  reverse_iterator rbegin() const {
    reverse_iterator it(head_);
    return it;
  }

  const_reverse_iterator crbegin() const {
    const_reverse_iterator it(head_);
    return it;
  }

  iterator end() const {
    iterator it(head_);
    return it;
  }

  const_iterator cend() const {
    const_iterator it(head_->next);
    return it;
  }

  reverse_iterator rend() const {
    reverse_iterator it(head_->next);
    return it;
  }

  const_reverse_iterator crend() const {
    const_reverse_iterator it(head_->next);
    return it;
  }

  void push_front(const T& value) {
    insert_(head_, value);
  }

  void push_back(const T& value) {
    insert_(head_->prev, value);
  }

  void pop_front() {
    erase_(head_->next);
  }

  void pop_back() {
    erase_(head_->prev);
  }

  void insert(const_iterator it, const T& value) {
    insert_(const_cast<Node*>((--it).it), value);
  }

  void erase(const_iterator it) {
    erase_(const_cast<Node*>(it.it));
  }
};
