# Projects on C++ language

* ### String

* ### Big Integer and Rational
  * BigInteger
    * fast multiplication using FFT
  * Rational

* ### Geometry

* ### Residue

* ### Matrix
  * Matrix
    * fast matrix multiplication using Strauss algorithm 

* ### Deque
  * Deque (simple version of std::deque)
    * random access iterator supporting
    * doesn't invalidate iterators

* ### List with Fast Allocator
  * FixedAllocator (singleton class allocating blocks of fixed size)
  * FastAllocator (stl compatible allocator)
    * 100% faster than std::allocator (in average)
  * List
    * allocators supporting (including allocators from std)
    * bidirectional iterators supporting
    * doesn't invalidate links and iterators

* ### Unordered Map
  * UnorderedMap (like std::unordered_map)
    * move semantics supporting
    * allocators supporting (including allocators from std)
    * forward iterators supporting
    * doesn't invalidate iterators and links
    * rehash O(n)

* ### Smart Pointers
  * SharedPtr (like std::shared_ptr)
    * Optimized memory allocations (implementation with Control Block)
    * move semantics supporting
    * allocators and deleters supporting (including allocators and deleters from std)
  * WeakPtr (like std::weak_ptr)
    * all features from SharedPtr
 
