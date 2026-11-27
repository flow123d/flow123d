#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <memory_resource>
#include <vector>
#include <iostream>
#include <new>
#include <stdexcept>   // !! Use Flow exception mechanism

#include "system/asserts.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"

#include <Eigen/Core>
#include <Eigen/Dense>


//TEST(FemToolsTest, working) {
//    size_t size = 256 * 1024;
//    std::shared_ptr<char> data_ptr(new char[size], std::default_delete<char[]>());
//    void* void_ptr = static_cast<void*>(data_ptr.get());
//
//    std::cout << "Test asm arena 1" << std::endl;
//    AssemblyArena asm_arena(void_ptr, size, 32);
//    std::cout << "Test asm arena 2" << std::endl;
//    uint data_size = 1212;
//
//    Eigen::Matrix<ArenaVec<double>, 3, 3> mat33_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 1> vec3_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 1> result_multi_arena;
//
//    for (uint i=0; i<9; ++i) {
//        mat33_arena(i) = ArenaVec<double>(data_size, asm_arena);
//        for (uint j=0; j<data_size; ++j) {
//            mat33_arena(i)(j) = 0.1 * (2*i + j + 1);
//        }
//    }
//    for (uint i=0; i<3; ++i) {
//        vec3_arena(i) = ArenaVec<double>(data_size, asm_arena);
//        for (uint j=0; j<data_size; ++j) {
//            vec3_arena(i)(j) = 0.1 * (2*i + j + 1);
//        }
//    }
//
//    std::cout << "Test patch arena 1" << std::endl;
//    PatchArena *patch_arena = asm_arena.get_child_arena();
//    std::cout << "Test patch arena 2" << std::endl;
//    for (uint i=0; i<9; ++i)
//        mat33_arena(i).set_patch_arena(*patch_arena);
//    for (uint i=0; i<3; ++i)
//        vec3_arena(i).set_patch_arena(*patch_arena);
//
//    for (uint i=0; i<1; ++i) {
//        patch_arena->reset();
//        result_multi_arena = mat33_arena * vec3_arena;
//    }
//}
//

/*
[100%] Linking CXX executable arena_test_bin
[100%] Built target arena_test_bin
[----------] Global test environment set-up.
[----------] 1 test from FemToolsTest
[ RUN      ] FemToolsTest.working
Test asm arena 1
Constructor 2: 262144, 32
Test asm arena 2
raw_aloocate used size: 16
raw_aloocate used size: 9712
raw_aloocate used size: 19408
raw_aloocate used size: 29104
raw_aloocate used size: 38800
raw_aloocate used size: 48496
raw_aloocate used size: 58192
raw_aloocate used size: 67888
raw_aloocate used size: 77584
raw_aloocate used size: 87280
raw_aloocate used size: 96976
raw_aloocate used size: 106672
Test patch arena 1
raw_aloocate used size: 116368
get_child_arena: 262144, 116368, 145776, 32
Constructor 2: 145776, 32
--------------------------------------------------------
Program Error: Violated assert!
> In file: /c/Users/davidflanderka/vyvoj/flow123d/src/fem/arena_resource.hh(70): Throw in function ArenaResource
> Expression: '(buffer_size%simd_alignment) == 0'
> Values:
  (buffer_size%simd_alignment) : '16'
  0 : '0'

** Stacktrace **

--------------------------------------------------------

===================================================================================
=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
=   PID 110 RUNNING AT 10cb26af67a8
=   EXIT CODE: 6
=   CLEANING UP REMAINING PROCESSES
=   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
===================================================================================
YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Aborted (signal 6)
This typically refers to a problem with your application.
Please see the FAQ page for debugging suggestions
make[3]: *** [unit_tests/system/CMakeFiles/arena-1-test.dir/build.make:70: unit_tests/system/CMakeFiles/arena-1-test] Error 6
make[2]: *** [CMakeFiles/Makefile2:7427: unit_tests/system/CMakeFiles/arena-1-test.dir/all] Error 2
make[1]: *** [CMakeFiles/Makefile2:7434: unit_tests/system/CMakeFiles/arena-1-test.dir/rule] Error 2
make: *** [Makefile:886: unit_tests/system/CMakeFiles/arena-1-test.dir/rule] Error 2
rel 4.0.3 /c/Users/davidflanderka/vyvoj/flow123d/unit_tests/system make -j 5 arena-1-test

POST-CHECKOUT(DF_patch_fe_data_tables -> DF_patch_fe_data_tables):
[  0%] Built target update_source_tree
[  0%] Built target gtest
[ 16%] Built target gtest_main
[ 83%] Built target system_lib
[ 83%] Built target gtest_mpi_obj
Consolidate compiler generated dependencies of target arena_test_bin
[100%] Building CXX object unit_tests/system/CMakeFiles/arena_test_bin.dir/arena_test.cpp.o
In file included from /c/Users/davidflanderka/vyvoj/flow123d/unit_tests/flow_gtest_mpi.hh:13,
                 from /c/Users/davidflanderka/vyvoj/flow123d/unit_tests/system/arena_test.cpp:3:
/c/Users/davidflanderka/vyvoj/flow123d/unit_tests/_gtest_extensions/gtest_throw_what.hh:22:22: warning: ‘std::string FormatDeathTestOutput(const string&)’ defined but not used [-Wunused-function]
   22 | static ::std::string FormatDeathTestOutput(const ::std::string& output) {
      |                      ^~~~~~~~~~~~~~~~~~~~~
[100%] Linking CXX executable arena_test_bin
[100%] Built target arena_test_bin
[----------] Global test environment set-up.
[----------] 1 test from FemToolsTest
[ RUN      ] FemToolsTest.working
Test asm arena 1
Constructor 1: 262144, 32
Test asm arena 2
raw_aloocate used size: 0
raw_aloocate used size: 9696
raw_aloocate used size: 19392
raw_aloocate used size: 29088
raw_aloocate used size: 38784
raw_aloocate used size: 48480
raw_aloocate used size: 58176
raw_aloocate used size: 67872
raw_aloocate used size: 77568
raw_aloocate used size: 87264
raw_aloocate used size: 96960
raw_aloocate used size: 106656
Test patch arena 1
raw_aloocate used size: 116352
get_child_arena: 262144, 116352, 145792, 32
Constructor 2: 145792, 32
Test patch arena 2
raw_aloocate used size: 0
raw_aloocate used size: 9696
raw_aloocate used size: 19392
raw_aloocate used size: 29088
raw_aloocate used size: 38784
raw_aloocate used size: 48480
raw_aloocate used size: 58176
raw_aloocate used size: 67872
raw_aloocate used size: 77568
raw_aloocate used size: 87264
raw_aloocate used size: 96960
raw_aloocate used size: 106656
raw_aloocate used size: 116352
raw_aloocate used size: 126048
raw_aloocate used size: 135744
[       OK ] FemToolsTest.working (0 ms)
[----------] 1 test from FemToolsTest (0 ms total)

[----------] Global test environment tear-down
[==========] 1 test from 1 test case ran. (0 ms total)
[  PASSED  ] 1 test.
[100%] Built target arena-1-test
 */

/// First proposal in issue
//class ArenaAllocatorTest : public std::pmr::memory_resource {
//public:
//    explicit ArenaAllocatorTest(size_t buffer_size)
//        : buffer(new char[buffer_size]),
//          buffer_size(buffer_size),
//          resource(buffer.get(), buffer_size, std::pmr::null_memory_resource()) {}
//          // !! use null_memory_resource in Debug to throw exception in case of out of arena space
//         // !! use  get_default_resource in Release (allow allocations out of arena to be from heap
//
//
//    ~ArenaAllocatorTest() = default;
//
//    // Set a checkpoint
//    void* checkpoint() {
//        return resource.upstream_resource();
//        //return resource.get_current_buffer();
//    }
//
//    // Reset to a checkpoint
//    void reset(void* checkpoint) {
//        resource.release();
//        //resource.set_current_buffer(checkpoint);
//    }
//
//protected:
//    // Override do_allocate to handle allocation logic
//    void* do_allocate(size_t bytes, size_t alignment) override {
//        void* p = resource.allocate(bytes, alignment);
//        if (p == nullptr) {  // test only in Debug when null_pointer_resource is in use
//            throw std::bad_alloc();
//        }
//        return p;
//    }
//
//    // Override do_deallocate (no-op for monotonic buffer)
//    void do_deallocate(void* p, size_t bytes, size_t alignment) override {
//        // No-op
//    }
//
//    // Override do_is_equal for memory resource comparison
//    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
//        return this == &other;
//    }
//
//private:
//    std::unique_ptr<char[]> buffer;
//    size_t buffer_size;
//    std::pmr::monotonic_buffer_resource resource;
//};
//
//TEST(Arena, chat_gpt_proposal_1) {
//    try {
//        ArenaAllocatorTest allocator(1024); // Create an arena with 1KB buffer
//
//        void* checkpoint1 = allocator.checkpoint(); // Set a checkpoint
//
//        void* p1 = allocator.allocate(256); // Allocate 256 bytes
//        void* p2 = allocator.allocate(512); // Allocate 512 bytes
//
//        allocator.reset(checkpoint1); // Reset to checkpoint, p1 and p2 are now invalid
//
//        // Allocate again after reset
//        void* p3 = allocator.allocate(128); // Allocate 128 bytes
//    } catch (const std::bad_alloc& e) {
//        // Handle allocation failure
//        std::cerr << "Allocation failed: " << e.what() << std::endl;
//    }
//
//}


/// Second proposal in issue
//class AlignedMemoryResourceTest : public std::pmr::memory_resource {
//public:
//    explicit AlignedMemoryResourceTest(std::pmr::monotonic_buffer_resource& upstream, size_t alignment)
//        : upstream_(upstream), alignment_(alignment) {}
//
//protected:
//    void* do_allocate(size_t bytes, size_t alignment) override {
//        //if (alignment_ > alignment) { // no effect, allocate methods aligns automatically
//        //    alignment = alignment_;
//        //}
//        return upstream_.allocate(bytes, alignment);
//    }
//
//    void do_deallocate(void* p, size_t bytes, size_t alignment) override {
//        upstream_.deallocate(p, bytes, alignment);
//    }
//
//    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
//        return this == &other;
//    }
//
//private:
//    std::pmr::monotonic_buffer_resource& upstream_;
//    size_t alignment_;
//};
//
//TEST(Arena, chat_gpt_proposal_2) {
//    // Shared buffer for the arena
//    const size_t buffer_size = 1024 * 1024; // 1MB buffer
//    char buffer[buffer_size];
//
//    // Create a monotonic buffer resource using the shared buffer
//    std::pmr::monotonic_buffer_resource shared_arena(buffer, buffer_size);
//
//    // Create aligned memory resources using the shared arena
//    AlignedMemoryResourceTest aligned_resource_8(shared_arena, 8);
//    AlignedMemoryResourceTest aligned_resource_256(shared_arena, 256);
//
//    // Create std::pmr::vector instances using the aligned memory resources
//    std::pmr::vector<int> vector_8(&aligned_resource_8);
//    std::pmr::vector<int> vector_256(&aligned_resource_256);
//    std::pmr::vector<double> vector_dbl_8(&aligned_resource_8);
//    std::pmr::vector<double> vector_dbl_256(&aligned_resource_256);
//
//    // Fill the vectors with some data
//    for (int i = 0; i < 10; ++i) {
//        vector_8.push_back(i);
//        vector_256.push_back(i * 10);
//        vector_dbl_8.push_back(i + 0.1);
//        vector_dbl_256.push_back(i * 10 + 0.5);
//    }
//
//    // Print the vectors
//    std::cout << "Vector int with 8-byte alignment: ";
//    for (int value : vector_8) {
//        std::cout << value << " ";
//    }
//    std::cout << std::endl;
//
//    std::cout << "Vector int with 256-byte alignment: ";
//    for (int value : vector_256) {
//        std::cout << value << " ";
//    }
//    std::cout << std::endl;
//
//    std::cout << "Vector double with 8-byte alignment: ";
//    for (double value : vector_dbl_8) {
//        std::cout << value << " ";
//    }
//    std::cout << std::endl;
//
//    std::cout << "Vector double with 256-byte alignment: ";
//    for (double value : vector_dbl_256) {
//        std::cout << value << " ";
//    }
//    std::cout << std::endl;
//}

namespace ArenaTest {

template <class Resource> class ArenaAllocator; // forward declaration


// Helper class of ArenaAllocator
template <class Resource>
class AlignedMemoryResource : public std::pmr::memory_resource {
public:
    explicit AlignedMemoryResource(std::pmr::monotonic_buffer_resource& upstream, size_t alignment)
        : upstream_(upstream), alignment_(alignment) {}

protected:
    void* do_allocate(size_t bytes, size_t alignment) override {
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
    Resource& upstream_;
    size_t alignment_;

    friend class ArenaAllocator<Resource>;
};


// Final proposal of Arena
template <class Resource>
class ArenaAllocator : public std::pmr::memory_resource {
public:
    explicit ArenaAllocator(size_t buffer_size, size_t simd_alignment)
        : buffer(new char[buffer_size]),
          buffer_size(buffer_size),
#ifdef FLOW123D_DEBUG
          resource_(buffer.get(), buffer_size, std::pmr::null_memory_resource()),
#else
          resource_(buffer.get(), buffer_size, std::pmr::get_default_resource()),
#endif
		  arena_resource_8_(this->resource_, 8),
		  arena_resource_simd_(this->resource_, simd_alignment),
		  simd_alignment_(simd_alignment)
    {}


    ~ArenaAllocator() = default;

    /// Getter for resource
    Resource &resource() {
    	return resource_;
    }


    /// Allocate and return data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    T* allocate_8(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_allocate with alignment argument
        return (T*)arena_resource_8_.do_allocate(bytes, arena_resource_8_.alignment_);
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    T* allocate_simd(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_allocate with alignment argument
        return (T*)arena_resource_simd_.do_allocate(bytes, arena_resource_simd_.alignment_);
    }

    /// Deallocate passed data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    void deallocate_8(T * data_ptr, size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_deallocate with alignment argument
        arena_resource_8_.do_deallocate((void *)data_ptr, bytes, arena_resource_8_.alignment_);
    }

    /// Deallocate passed data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    void deallocate_simd(T * data_ptr, size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_deallocate with alignment argument
        arena_resource_simd_.do_deallocate((void *)data_ptr, bytes, arena_resource_simd_.alignment_);
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
    AlignedMemoryResource<Resource> arena_resource_8_;
    AlignedMemoryResource<Resource> arena_resource_simd_;
    size_t simd_alignment_;
};


template<class T>
class ArenaVec {
public:
    /// Type definitions
	typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayData;
    typedef Eigen::Vector<ArrayData, Eigen::Dynamic> TableData;
    typedef Eigen::Matrix<ArrayData, Eigen::Dynamic, Eigen::Dynamic> MatrixData;

    /// Base constructor, set invalid data pointer
    ArenaVec(std::vector<size_t> shape)
    : shape_( set_shape_vec(shape) ), row_size_(0), data_ptr_(nullptr) {}

    /**
     * Constructor. Set sizes and allocate data pointer
     */
    ArenaVec(std::vector<size_t> shape, size_t row_size, ArenaAllocator<std::pmr::monotonic_buffer_resource> &arena)
    : shape_( set_shape_vec(shape) ), data_ptr_(nullptr) {
        allocate_data(row_size, arena);
    }

    /// Return number of components
    inline size_t n_comp() const {
        return shape_[0] * shape_[1];
    }

    /**
     * Allocate data on \p arena memory
     *
     * Size of data is given by \p shape_ and \p row_size
     */
    void allocate_data(size_t row_size, ArenaAllocator<std::pmr::monotonic_buffer_resource> &arena) {
        if (data_ptr_ != nullptr) {
            arena.deallocate_simd<T>(data_ptr_, n_comp() * row_size );
        }
        row_size_ = row_size;
        data_ptr_ = arena.allocate_simd<T>( n_comp() * row_size );
    }

    /**
     * Maps data pointer to Eigen Map of dimensions given by shape_ and row_size_ and returns it.
     */
    Eigen::Map<MatrixData> data_map() {
        // Create the TableDbl
        TableData table_data( n_comp() );

        // Map each segment of the array d to the inner ArrayDbl of TableDbl
        for (int i = 0; i < n_comp(); ++i) {
            table_data(i) = Eigen::Map<ArrayData>(data_ptr_ + i * row_size_, row_size_);
        }

        return Eigen::Map<MatrixData>(table_data.data(), shape_[0], shape_[1]);
    }

//    Eigen::Map<MatrixData> full_map() {
//        Eigen::Map< Eigen::Matrix<ArrayData, Eigen::Dynamic, Eigen::Dynamic> > map(data_ptr_, shape_[0], shape_[1]);
//        return map;
//    }

    /// Temporary test
    MatrixData data_matrix() {
        // Create the TableDbl
    	MatrixData matrix_data( shape_[0], shape_[1] );

        // Map each segment of the array d to the inner ArrayDbl of TableDbl
        for (int i = 0; i < n_comp(); ++i) {
        	matrix_data(i%shape_[0], i/shape_[0]) = Eigen::Map<ArrayData>(data_ptr_ + i * row_size_, row_size_);
        }

        return matrix_data;
    }

    /// Return data pointer (development method)
    T* data_ptr() {
        return data_ptr_;
    }

protected:
    std::vector<size_t> set_shape_vec(std::vector<size_t> shape) const {
    	ASSERT( (shape.size() == 1) || (shape.size() == 2) )(shape.size()).error("Invalid size of shape vector!\n");
    	if (shape.size() == 1) shape.push_back(1);
    	return shape;
    }

    std::vector<size_t> shape_;
    size_t row_size_;
    T* data_ptr_;
};

}; // namespace ArenaTest

TEST(Arena, allocatzor) {
    try {
        ArenaTest::ArenaAllocator<std::pmr::monotonic_buffer_resource> allocator(1024 * 1024, 256); // Create an arena with 1MB buffer

        // Create std::pmr::vector instances using the aligned memory resources
        ArenaTest::AlignedMemoryResource<std::pmr::monotonic_buffer_resource> aligned_resource_8(allocator.resource(), 8);
        ArenaTest::AlignedMemoryResource<std::pmr::monotonic_buffer_resource> aligned_resource_simd(allocator.resource(), 256);
        std::pmr::vector<int> vector_8( &aligned_resource_8 );
        std::pmr::vector<int> vector_256( &aligned_resource_simd );
        std::pmr::vector<double> vector_dbl_8( &aligned_resource_8 );
        std::pmr::vector<double> vector_dbl_256( &aligned_resource_simd );
        // Create array of double of length 10
        double *custom_memory = allocator.allocate_simd<double>(10);
        ArenaTest::ArenaVec<double> quantity_deta({2, 3}, 10, allocator);

        // Fill the vectors with some data
        for (int i = 0; i < 10; ++i) {
            vector_8.push_back(i);
            vector_256.push_back(i * 10);
            vector_dbl_8.push_back(i + 0.1);
            vector_dbl_256.push_back(i * 10 + 0.5);
            custom_memory[i] = i * 10 + 0.5;
        }
        double* data_ptr = quantity_deta.data_ptr();
        for (int i = 0; i < 60; ++i) {
            data_ptr[i] = (i%10 + 1) + 0.1 * (i/10 + 1);
        }

        // Print the vectors
        std::cout << "Vector int with 8-byte alignment: ";
        for (int value : vector_8) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::cout << "Vector int with 256-byte alignment: ";
        for (int value : vector_256) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::cout << "Vector double with 8-byte alignment: ";
        for (double value : vector_dbl_8) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::cout << "Vector double with 256-byte alignment: ";
        for (double value : vector_dbl_256) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::cout << "Array of double with 256-byte alignment: ";
        for (uint i=0; i<10; ++i) {
            std::cout << custom_memory[i] << " ";
        }
        std::cout << std::endl << std::endl;

        // Map the custom allocated memory to an Eigen vector
        Eigen::Map<Eigen::VectorXd> eigen_vector(custom_memory, 10);
        Eigen::Map<Eigen::MatrixXd> eigen_matrix(custom_memory, 2, 5);
        // ... or Eigen::Map<Eigen::MatrixXd> eigen_matrix(&vector_dbl_256[0], 2, 5);

        // Print the Eigen vector and matrix
        std::cout << "Eigen vector:\n" << eigen_vector.transpose() << std::endl;
        std::cout << "Eigen matrix:\n" << eigen_matrix << std::endl;

        // Change value in vector, print matrix, value (1,0) must be changed
        eigen_vector(1) = 5.5;
        std::cout << "Eigen matrix changed:\n" << eigen_matrix << std::endl << std::endl;

//        auto map = quantity_deta.data_map();
//        std::cout << "ArenaVec dimensions: " << map.rows() << ", " << map.cols() << std::endl;

        auto matrix = quantity_deta.data_matrix();
        std::cout << "ArenaVec dimensions: " << matrix.rows() << ", " << matrix.cols() << ", " << matrix(0,0).rows() << std::endl;
        std::cout << "Matrix 2x3:" << std::endl;
        for (int i = 0; i < 6; ++i) {
        	std::cout << "(" << i%2 << "," << i/2 << "): ";
        	std::cout << matrix(i%2, i/2).transpose() << std::endl;
        }
        quantity_deta.data_ptr()[1] = 0.5;
        std::cout << "Matrix 2x3 changed:" << std::endl;
        for (int i = 0; i < 6; ++i) {
        	std::cout << "(" << i%2 << "," << i/2 << "): ";
        	std::cout << matrix(i%2, i/2).transpose() << std::endl;
        }
    } catch (const std::bad_alloc& e) {
        // Handle allocation failure
        std::cerr << "Allocation failed: " << e.what() << std::endl;
    }

}

TEST(Arena, map_access_speed_test) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    static const uint n_loops = 1e06;
    size_t single_sum=0, repeated_sum=0;

    ArenaTest::ArenaAllocator<std::pmr::monotonic_buffer_resource> allocator(1024 * 1024, 256); // Create an arena with 1MB buffer

    // Create std::pmr::vector instances using the aligned memory resources
    ArenaTest::ArenaVec<double> quantity_deta({2, 3}, 10, allocator);

    // Fill the vectors with some data
    double* data_ptr = quantity_deta.data_ptr();
    for (int i = 0; i < 60; ++i) {
        data_ptr[i] = (i%10 + 1) + 0.1 * (i/10 + 1);
    }

    START_TIMER("single_map_construct");
    auto map_1 = quantity_deta.data_map();
    for (uint i=0; i<n_loops; ++i) {
        if (i%2==0) single_sum += map_1.rows();
        else single_sum += map_1.cols();
    }
    END_TIMER("single_map_construct");

    START_TIMER("repeated_map_construct");
    for (uint i=0; i<n_loops; ++i) {
        auto map_2 = quantity_deta.data_map();
        if (i%2==0) repeated_sum += map_1.rows();
        else repeated_sum += map_1.cols();
    }
    END_TIMER("repeated_map_construct");
    EXPECT_EQ(single_sum, repeated_sum);

    // profiler autput
    FilePath fp("map_access_profiler.json", FilePath::output_file);
	Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
}

/*
 * Result:
 *
 * n_repeats                     1E05        1E06
 * single_map_construct     0.0001052   0.0006767
 * repeated_map_construct   0.0291635   0.2823245
 */


TEST(Arena, asm_patch_model) {
	AssemblyArena asm_arena(1024 * 1024, 256);

//    ArenaVec<double> data_vec(10, asm_arena);
//    double *custom_memory1 = data_vec.data_ptr();
//
//    for (int i = 0; i < 10; ++i) {
//        custom_memory1[i] = i * 10 + 0.5;
//    }
//    std::cout << "data_vec: " << data_vec.eigen_map().transpose() << std::endl;

    Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac(2,3);
    for (uint i=0; i<2*3; ++i)
        jac(i) = ArenaVec<double>(10, asm_arena);
    std::cout << "jac: " << jac.rows() << " - " << jac.cols() << std::endl;
    std::cout << "inner: " << jac(1,2).data_size() << std::endl;
    jac(0,0)(0) = 1.5;
    std::cout << "zero item: " << jac(0,0)(0) << std::endl;
    Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > submat(jac.data(), 2, 2);
    std::cout << "mat zero item: " << submat(0,0)(0) << std::endl;
    jac(0,0)(0) = 2.5;
    std::cout << "changed zero item: " << jac(0,0)(0) << std::endl;
    std::cout << "changed mat zero item: " << submat(0,0)(0) << std::endl;

}


TEST(Arena, arena_alloc) {
    static const uint N_RUNS = 4e5;

    AssemblyArena asm_arena(64 * 1024, 32, std::pmr::get_default_resource());
    size_t sum_sizes = 0;

    ArenaVec<double> data_vec1(4 * 1024 + 4, asm_arena);
    sum_sizes += data_vec1.data_size();

    PatchArena *patch_arena = asm_arena.get_child_arena();

    for (uint i=0; i<N_RUNS; ++i) {
        patch_arena->reset();
        ArenaVec<double> data_vec2(4 * 1024 + 4, *patch_arena);
        sum_sizes += data_vec2.data_size();
    }
    std::cout << "sum_sizes: " << sum_sizes << std::endl;
}


/**
 * Tests deallocation of arena memory.
 *
 * Allocate arena and 3 vectors, first vector is small and is stored in AsmArena.
 * Two remaining vectors are large and are stored in PatchArena. Size of these vector is set up like this:
 * 1. All 3 vectors are stored in resource of arena.
 * 2. Sum of sizes of all vectors is larger than size of arena and allocation needs to get second block of memory.
 *    If deallocation is not call, test causes memory error.
 */
TEST(Arena, dealloc_test) {
    static const uint N_RUNS = 4e6;
    static const size_t arena_vec_size = 1024;
    static const size_t patch_vec_size = 510 * 1024; // first case - see test description
    //static const size_t patch_vec_size = 512 * 1024; // second case

    AssemblyArena asm_arena(1024 * 1024, 256, std::pmr::get_default_resource());
    size_t sum_sizes = 0;

    ArenaVec<double> data_vec1(arena_vec_size, asm_arena);
    sum_sizes += data_vec1.data_size();

    PatchArena *patch_arena = asm_arena.get_child_arena();

    for (uint i=0; i<N_RUNS; ++i) {
        patch_arena->reset();
        ArenaVec<double> data_vec2(patch_vec_size, *patch_arena);
        ArenaVec<double> data_vec3(patch_vec_size, *patch_arena);
        sum_sizes += data_vec2.data_size() + data_vec3.data_size();
    }
    std::cout << "sum_sizes: " << sum_sizes << std::endl;
}

