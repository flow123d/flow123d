#ifndef ARENA_RESOURCE_HH_
#define ARENA_RESOURCE_HH_

#include <memory_resource>
#include <vector>
#include <iostream>
#include <new>
#include <stdexcept>   // !! Use Flow exception mechanism

#include "system/asserts.hh"

#include <Eigen/Core>
#include <Eigen/Dense>


template <class Resource> class ArenaResource; // forward declaration


// Helper class of ArenaResource, allows allocate aligned blocks of data
template <class Upstream>
class AlignedMemoryResource : public std::pmr::memory_resource {
public:
    explicit AlignedMemoryResource(std::pmr::monotonic_buffer_resource& upstream, size_t alignment)
        : upstream_(upstream), alignment_(alignment) {}

protected:
    void* do_allocate(size_t bytes, size_t alignment) override {
        if (alignment_ > alignment) {
            alignment = alignment_;
        }
        void* p = upstream_.allocate(bytes, alignment);
        if (p == nullptr) {  // test only in Debug when null_pointer_resource is in use
            throw std::bad_alloc();
        }
        return p;
    }

    void do_deallocate(void* p, size_t bytes, size_t alignment) override {
        upstream_.deallocate(p, bytes, alignment);
    }

    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
        return this == &other;
    }

private:
    Upstream& upstream_;
    size_t alignment_;

    friend class ArenaResource<Upstream>;
};


// Final proposal of Arena
template <class Resource>
class ArenaResource : public std::pmr::memory_resource {
public:
    explicit ArenaResource(size_t buffer_size, size_t simd_alignment)
        : buffer(new char[buffer_size]),
          buffer_size(buffer_size),
#ifdef FLOW123D_DEBUG
          resource_(buffer.get(), buffer_size, std::pmr::null_memory_resource()),
#else
          resource_(buffer.get(), buffer_size, std::pmr::get_default_resource()),
#endif
          simd_alignment_(simd_alignment)
    {}


    ~ArenaResource() = default;

    /// Getter for resource
    Resource &resource() {
    	return resource_;
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    T* allocate_8(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        AlignedMemoryResource<Resource> a_res(this->resource_, 8);
        return (T*)a_res.do_allocate(bytes, a_res.alignment_);
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    T* allocate_simd(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        AlignedMemoryResource<Resource> a_res(this->resource_, simd_alignment_);
        return (T*)a_res.do_allocate(bytes, a_res.alignment_);
    }

    // Reset allocated data
    void reset() {
        resource_.release();
    }

protected:
    /// Override do_allocate to handle allocation logic
    void* do_allocate(size_t bytes, size_t alignment) override {
        // No-op
        return nullptr;
    }

    /// Override do_deallocate (no-op for monotonic buffer)
    void do_deallocate(void* p, size_t bytes, size_t alignment) override {
        // No-op
    }

    /// Override do_is_equal for memory resource comparison
    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
        return this == &other;
    }

private:
    std::unique_ptr<char[]> buffer;
    size_t buffer_size;
    Resource resource_;
    size_t simd_alignment_;
};

using AssemblyArena = ArenaResource<std::pmr::monotonic_buffer_resource>;
using PatchArenaResource = ArenaResource<AssemblyArena>;



#endif /* ARENA_RESOURCE_HH_ */
