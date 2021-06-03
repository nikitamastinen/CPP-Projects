//
// Created by Nikita Mastinen on 02.05.2021.
//

#include <type_traits>
#include <memory>

class ControlBlockBase {
public:
  size_t shared_counter = 0;
  size_t weak_counter = 0;

  explicit ControlBlockBase(size_t shared_counter) : shared_counter(shared_counter), weak_counter(0) {}

  virtual void destroy() = 0;

  virtual void destroy_value() = 0;

  virtual void destroy_value(void*) = 0;

protected:
  virtual ~ControlBlockBase() = default;
};

template<typename Y, typename Allocator, typename Deleter>
class EmptyControlBlock: public ControlBlockBase {
public:
  Allocator allocator;
  Deleter object_deleter;

  using SelfAllocator = typename std::allocator_traits<Allocator>
  ::template rebind_alloc<EmptyControlBlock>;

  explicit EmptyControlBlock(size_t shared_counter, const Allocator& allocator, const Deleter& y_deleter) :
      ControlBlockBase(shared_counter),
      allocator(allocator),
      object_deleter(y_deleter) {
  }

  void destroy_value() override {};

  void destroy_value(void* y) override {
    object_deleter.operator()(static_cast<Y*>(y));
  }

  void destroy() override {
    SelfAllocator copy_allocator(allocator);
    std::allocator_traits<SelfAllocator>::deallocate(copy_allocator, this, 1);
  }

  ~EmptyControlBlock() override {}
};

template<
    typename T,
    typename Allocator = std::allocator<T>
>
class ControlBlock: public ControlBlockBase {
public:
  T value;
  Allocator object_allocator;
  using SelfAllocator = typename std::allocator_traits<Allocator>
  ::template rebind_alloc<ControlBlock>;

  template<typename ...Args>
  explicit ControlBlock(const Allocator& object_allocator, Args&&... args)
      : ControlBlockBase(1),
        value(std::forward<Args>(args)...),
        object_allocator(object_allocator)
  {}

  void destroy_value() override {
    std::allocator_traits<Allocator>::destroy(object_allocator, &value);
  }

  void destroy_value(void*) override {}

  void destroy() override {
    SelfAllocator copy_allocator(object_allocator);
    object_allocator.~Allocator();
    std::allocator_traits<SelfAllocator>::deallocate(copy_allocator, this, 1);
  }

  ~ControlBlock() override {}
};

template<typename T>
class WeakPtr;

template<typename T>
class SharedPtr {
  template<typename Y>
  friend class SharedPtr;

  template<typename Y>
  friend class WeakPtr;

  template <typename Y, typename Allocator, typename... Args>
  friend SharedPtr<Y> allocateShared(Allocator, Args&&...);
private:
  T* value;
  ControlBlockBase* counter;

  //-----------------------------------CONSTRUCTORS------------------------------------------------
  template<typename Allocator>
  explicit SharedPtr(ControlBlock<T, Allocator>* other, int)
      : value(nullptr),
        counter(other)
  {}
public:
  SharedPtr() : value(nullptr), counter(nullptr) {}

  template<typename Y, typename Deleter, typename Allocator>
  SharedPtr(Y* object, const Deleter& deleter, const Allocator& allocator) : value(object)  {
    using ControlBlockType = EmptyControlBlock<Y, Allocator, Deleter>;
    using BlockAllocatorType =
    typename std::allocator_traits<Allocator>::template rebind_alloc<ControlBlockType>;
    BlockAllocatorType block_allocator(allocator);
    ControlBlockType* new_counter =
        std::allocator_traits<BlockAllocatorType>::allocate(block_allocator, 1);
    ::new (new_counter) ControlBlockType(1, allocator, deleter);
    counter = new_counter;
  }

  template<typename Y, typename Deleter>
  SharedPtr(Y* object, const Deleter& deleter) : SharedPtr(object, deleter, std::allocator<Y>()) {}

  template<typename Y>
  SharedPtr(Y* object) : SharedPtr(object, std::default_delete<Y>(), std::allocator<Y>()) {}

  SharedPtr(const SharedPtr& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->shared_counter;
  }

  template<typename Y>
  SharedPtr(const SharedPtr<Y>& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->shared_counter;
  }

  template<typename Y>
  SharedPtr(SharedPtr<Y>&& other) :
      value(std::move(other.value)),
      counter(std::move(other.counter)) {
    other.counter = nullptr;
    other.value = nullptr;
  }

  template <typename Y>
  explicit SharedPtr(const WeakPtr<Y>& other) : value(other.value), counter(other.counter) {
    ++counter->shared_counter;
  }

  //-----------------------------------OPERATORS--------------------------------------------------
  template<typename Y>
  SharedPtr& operator=(const SharedPtr<Y>& other) {
    SharedPtr copy = other;
    swap(copy);
    return *this;
  }

  SharedPtr& operator=(const SharedPtr<T>& other) {
    if (&other == this) {
      return *this;
    }
    SharedPtr copy = other;
    swap(copy);
    return *this;
  }

  SharedPtr& operator=(SharedPtr<T>&& other) noexcept {
    SharedPtr moved = std::move(other);
    swap(moved);
    return *this;
  }

  template<typename Y>
  SharedPtr& operator=(SharedPtr<Y>&& other) noexcept {
    SharedPtr moved = std::move(other);
    swap(moved);
    return *this;
  }

  T& operator*() const {
    if (value == nullptr) {
      return static_cast<ControlBlock<T>*>(counter)->value;
    }
    return *value;
  }

  T* operator->() const {
    return &(**this);
  }

  T* get() const {
    return value;
  }

  //----------------------------------------------------------------------------------------------

  size_t use_count() const {
    return counter->shared_counter;
  }

  void swap(SharedPtr<T>& other) {
    std::swap(value, other.value);
    std::swap(counter, other.counter);
  }

  template<typename Y>
  void reset(Y* y) {
    SharedPtr(y).swap(*this);
  }

  void reset() {
    SharedPtr().swap(*this);
  }

  //------------------------------------------------------------------------------------------------

  ~SharedPtr() {
    if (counter == nullptr) return;
    --counter->shared_counter;
    if (counter->shared_counter != 0) {
      return;
    }
    if (value == nullptr) {
      counter->destroy_value();
    } else {
      counter->destroy_value(value);
    }
    if (counter->weak_counter == 0) counter->destroy();
  }
};

template <typename T, typename Allocator, typename... Args>
SharedPtr<T> allocateShared(Allocator allocator, Args&&... args) {
  using ControlBlockType = ControlBlock<T, Allocator>;
  using BlockAllocatorType =
  typename std::allocator_traits<Allocator>::template rebind_alloc<ControlBlockType>;
  BlockAllocatorType block_allocator;
  ControlBlockType* ptr =
      std::allocator_traits<BlockAllocatorType>::allocate(block_allocator, 1);
  std::allocator_traits<BlockAllocatorType>::construct(
      block_allocator, ptr, allocator, std::forward<Args>(args)...
  );
  return SharedPtr<T>(ptr, 1);
}


template <typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
  return allocateShared<T>(std::allocator<T>(), std::forward<Args>(args)...);
}

template <typename T>
class WeakPtr {
  template<typename Y>
  friend class WeakPtr;

  template<typename Y>
  friend class SharedPtr;
private:
  T* value;
  ControlBlockBase* counter;
public:
  //-----------------------------------CONSTRUCTORS------------------------------------------------
  WeakPtr() : value(nullptr), counter(nullptr) {}

  bool expired() const {
    return counter == nullptr || counter->shared_counter == 0;
  }

  SharedPtr<T> lock() const {
    return expired() ? SharedPtr<T>() : SharedPtr<T>(*this);
  }

  WeakPtr(const WeakPtr& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->weak_counter;
  }

  WeakPtr(const SharedPtr<T>& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->weak_counter;
  }

  template<typename Y>
  WeakPtr(const WeakPtr<Y>& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->weak_counter;
  }

  template<typename Y>
  WeakPtr(const SharedPtr<Y>& other) : value(other.value), counter(other.counter) {
    if (counter != nullptr) ++counter->weak_counter;
  }

  template<typename Y>
  WeakPtr(WeakPtr<Y>&& other)
      : value(std::move(other.value)),
        counter(std::move(other.counter)) {
    other.counter = nullptr;
    other.value = nullptr;
  }

  template<typename Y>
  WeakPtr(SharedPtr<Y>&& other)
      : value(std::move(other.value)),
        counter(std::move(other.counter)) {
    --counter->shared_counter;
    ++counter->weak_counter;
    other.counter = nullptr;
    other.value = nullptr;
  }

  //-----------------------------------OPERATORS--------------------------------------------------
  WeakPtr& operator=(const SharedPtr<T>& other) {
    WeakPtr copy = other;
    swap(copy);
    return *this;
  }

  template<typename Y>
  WeakPtr& operator=(const SharedPtr<Y>& other) {
    WeakPtr copy = other;
    swap(copy);
    return *this;
  }

  WeakPtr& operator=(SharedPtr<T>&& other) noexcept {
    WeakPtr moved = std::move(other);
    swap(moved);
    other.counter = nullptr;
    other.value = nullptr;
    return *this;
  }

  template<typename Y>
  WeakPtr& operator=(SharedPtr<Y>&& other) noexcept {
    WeakPtr moved = std::move(other);
    swap(moved);
    return *this;
  }
//-------------------------------------------------------------------------------------------------
  void swap(WeakPtr<T>& other) {
    std::swap(value, other.value);
    std::swap(counter, other.counter);
  }

  T& operator*() const {
    return *value;
  }

  T* operator->() const {
    return &(**this);
  }

  size_t use_count() const {
    return counter == nullptr ? 0 : counter->shared_counter;
  }

  ~WeakPtr() {
    if (counter == nullptr) return;
    --counter->weak_counter;
    if (counter->weak_counter > 0) {
      return;
    }
    if (counter->shared_counter == 0) {
      counter->destroy();
    }
  }
};