#include <iostream>
#include <chrono>
#include <functional>
#include "helpers.hpp"
#include "implementation.hpp"

int main(int argc, char **argv)
{
    bool print_runtime_only = false;
    int exec_mode = 0; // 0. all, 1 sequential only, 2. parallel only
	int grain_size = 256; // optional parameter to use for adjusting task granularity
	std::string X_filename = "data/X.txt", Y_filename = "data/Y.txt", output_filename = "data/aligned-sequential.txt";
	std::string substitution_matrix_file = "data/blosum62.txt";
	parse_args(argc, argv, X_filename, Y_filename, output_filename, grain_size, exec_mode, print_runtime_only);
    unsigned long entries_visited = 0, entries_visited_sequential = 0;

    SequenceInfo sinfo(X_filename, Y_filename);
    std::cout << "Loaded " << X_filename << " and " << Y_filename << " sequences with sizes " << sinfo.rows -1  << " and " << sinfo.cols -1 << std::endl;
    std::cout << "Matrix S size: [" << sinfo.rows << "x" << sinfo.cols << "]" << std::endl;

    // allocate
    float** S = allocate(sinfo.rows,sinfo.cols, 0); // Similarity Matrix

    std::unordered_map<char, int> cmap; // map Amino Acid (a character) to an index in Substitution Matrix

    float** SUB = sinfo.substitution_matrix_from_file(substitution_matrix_file, cmap);

    unsigned long expected_visited = (unsigned long)(sinfo.rows-1)*(sinfo.cols-1)+sinfo.rows+sinfo.cols-1;
    
    // sequential version
    if ( exec_mode == 1 || exec_mode < 1) {
        auto t_seq_1 = std::chrono::high_resolution_clock::now();

        entries_visited_sequential = sinfo.gpsa_sequential(S, SUB, cmap);
        
        auto t_seq_2 = std::chrono::high_resolution_clock::now();
        
        sinfo.traceback_and_save(output_filename, S, SUB, cmap, false);
        std::cout << "\n== Sequential version completed in " << std::chrono::duration<float>(t_seq_2 - t_seq_1).count() << " seconds." << std::endl; 
        std::cout << "   Entries visited: " << entries_visited_sequential << " " << (expected_visited ? "" : "NOT OK") << std::endl; 
        std::cout << "   Score: " << S[sinfo.rows-1][sinfo.cols-1] << ", Similarity Score: " << sinfo.similarity_score << ", Identity Score: " << sinfo.identity_score << ", Gaps: " << sinfo.gap_count << ", Length (with gaps): " << sinfo.X_aligned.size() << std::endl; 
        sinfo.reset(S);
    }
    
    // parallel version
    if ( exec_mode == 2 || exec_mode < 1) {
        auto t_taskloop_1 = std::chrono::high_resolution_clock::now();

        entries_visited = sinfo.gpsa_parallel(S, SUB, cmap, grain_size);
        
        auto t_taskloop_2 = std::chrono::high_resolution_clock::now();
        
        sinfo.traceback_and_save("data/aligned-taskloop.txt", S, SUB, cmap);

        std::cout << "\n== Parallel version completed in " << std::chrono::duration<float>(t_taskloop_2 - t_taskloop_1).count() << " seconds." << std::endl; 
        std::cout << "   Entries visited: " << entries_visited << " " << (expected_visited == entries_visited ? "" : "NOT OK") << std::endl; 
        std::cout << "   Score: " << S[sinfo.rows-1][sinfo.cols-1] << ", Similarity Score: " << sinfo.similarity_score << ", Identity Score: " << sinfo.identity_score << ", Gaps: " << sinfo.gap_count << ", Length (with gaps): " << sinfo.X_aligned.size() << std::endl; 
        if ( exec_mode < 1 ) std::cout << "   Checking results: " << (sinfo.verify(output_filename, "data/aligned-taskloop.txt") ? "OK" : "NOT OK") << std::endl;
        sinfo.reset(S);
    }

    deallocate(S);
    deallocate(SUB);

    return 0;
}
