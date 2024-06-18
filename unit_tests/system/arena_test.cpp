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

#include <Eigen/Core>
#include <Eigen/Dense>


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


template <class Resource> class ArenaAllocator; // forward declaration


// Helper class of ArenaAllocator
template <class Resource>
class AlignedMemoryResource : public std::pmr::memory_resource {
public:
    explicit AlignedMemoryResource(std::pmr::monotonic_buffer_resource& upstream, size_t alignment)
        : upstream_(upstream), alignment_(alignment) {}

protected:
    inline size_t align_size(size_t size) {
    	return (size + alignment_ - 1) / alignment_ * alignment_;
    }

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
          resource(buffer.get(), buffer_size, std::pmr::null_memory_resource()),
#else
          resource(buffer.get(), buffer_size, std::pmr::get_default_resource()),
#endif
		  aligned_resource_8_(this->resource, 8),
		  aligned_resource_simd_(this->resource, simd_alignment)
    {}


    ~ArenaAllocator() = default;

    /// Getter for aligned_resource_8_
    AlignedMemoryResource<Resource> aligned_resource_8() {
    	return aligned_resource_8_;
    }

    /// Getter for aligned_resource_simd_
    AlignedMemoryResource<Resource> aligned_resource_simd() {
    	return aligned_resource_simd_;
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    T* allocate_8(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_allocate with alignment argument
        return (T*)aligned_resource_8_.do_allocate(bytes, aligned_resource_8_.alignment_);
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    T* allocate_simd(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_allocate with alignment argument
        return (T*)aligned_resource_simd_.do_allocate(bytes, aligned_resource_simd_.alignment_);
    }

    /// Deallocate passed data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    void deallocate_8(T * data_ptr, size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_deallocate with alignment argument
        aligned_resource_8_.do_deallocate((void *)data_ptr, bytes, aligned_resource_8_.alignment_);
    }

    /// Deallocate passed data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    void deallocate_simd(T * data_ptr, size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        // Call override method do_deallocate with alignment argument
        aligned_resource_simd_.do_deallocate((void *)data_ptr, bytes, aligned_resource_simd_.alignment_);
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
    Resource resource;
    AlignedMemoryResource<Resource> aligned_resource_8_;
    AlignedMemoryResource<Resource> aligned_resource_simd_;
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


TEST(Arena, allocatzor) {
    try {
        ArenaAllocator<std::pmr::monotonic_buffer_resource> allocator(1024 * 1024, 256); // Create an arena with 1MB buffer

        // Create std::pmr::vector instances using the aligned memory resources
        AlignedMemoryResource<std::pmr::monotonic_buffer_resource> aligned_resource_8 = allocator.aligned_resource_8();
        AlignedMemoryResource<std::pmr::monotonic_buffer_resource> aligned_resource_simd = allocator.aligned_resource_simd();
        std::pmr::vector<int> vector_8( &aligned_resource_8 );
        std::pmr::vector<int> vector_256( &aligned_resource_simd );
        std::pmr::vector<double> vector_dbl_8( &aligned_resource_8 );
        std::pmr::vector<double> vector_dbl_256( &aligned_resource_simd );
        // Create array of double of length 10
        double *custom_memory = allocator.allocate_simd<double>(10);
        ArenaVec<double> quantity_deta({2, 3}, 10, allocator);

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

    ArenaAllocator<std::pmr::monotonic_buffer_resource> allocator(1024 * 1024, 256); // Create an arena with 1MB buffer

    // Create std::pmr::vector instances using the aligned memory resources
    ArenaVec<double> quantity_deta({2, 3}, 10, allocator);

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
