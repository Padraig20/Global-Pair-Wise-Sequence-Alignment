#include <unordered_map>
#include <omp.h>
#include "helpers.hpp"
#include <algorithm>

unsigned long SequenceInfo::gpsa_sequential(float** S, float** SUB, std::unordered_map<char, int>& cmap) {
    unsigned long visited = 0;
    gap_penalty = SUB[0][cmap['*']]; // min score

    unsigned long hey = 0;

	// Boundary
    for (int i = 1; i < rows; i++) {
        S[i][0] = i * gap_penalty;
		hey++;
	}

    for (int j = 0; j < cols; j++) {
        S[0][j] = j * gap_penalty;
		hey++;
	}

	// Main part
	for (int i = 1; i < rows; i++) {
		for (int j = 1; j < cols; j++) {
			float match = S[i - 1][j - 1] + SUB[ cmap.at(X[i - 1]) ][ cmap.at(Y[j-1]) ];
			float del = S[i - 1][j] + gap_penalty;
			float insert = S[i][j - 1] + gap_penalty;
			S[i][j] = std::max({match, del, insert});

			visited++;
		}
	}

	std::cout << hey << "\n";
	std::cout << visited << "\n";

    return visited + hey;
}

unsigned long SequenceInfo::gpsa_parallel(float** S, float** SUB, std::unordered_map<char, int> cmap, int grainsize = 1) {
    unsigned long visited = 0;
    gap_penalty = SUB[0][cmap['*']]; // min score

    std::vector<unsigned long> thread_counters(omp_get_max_threads(), 0);

	#pragma omp parallel
    {
		#pragma omp for
		for (int i = 1; i < rows; i++) {
			S[i][0] = i * gap_penalty;
			thread_counters[omp_get_thread_num()]++;
		}

		#pragma omp for
		for (int j = 0; j < cols; j++) {
			S[0][j] = j * gap_penalty;
			thread_counters[omp_get_thread_num()]++;
		}
    }

	#pragma omp parallel
    {
        #pragma omp single
        {
            int waves = rows + cols - 1;
            for (int wave = 0; wave < waves; wave++) {
            #pragma omp taskloop grainsize(grainsize)
                for (int i = 1; i < rows; i++) {
                    int j = wave - i;
                    if (j >= 1 && j < cols) {
						float match = S[i - 1][j - 1] + SUB[ cmap.at(X[i - 1]) ][ cmap.at(Y[j-1]) ];
						float del = S[i - 1][j] + gap_penalty;
						float insert = S[i][j - 1] + gap_penalty;
						S[i][j] = std::max({match, del, insert});
						#pragma omp atomic update
                        thread_counters[omp_get_thread_num()]++;
                    }
                }
            }
        }
    }

	std::cout << "Workload partitioning: ";

    for (unsigned long counter : thread_counters) {
            std::cout << counter << ' ';
    }

    std::cout << std::endl;

    for (unsigned long counter : thread_counters) {
        visited += counter;
	}

	std::cout << "Is:        " << visited << "\n";

    std::cout << "Should be: 350418241\n";

    return visited;
}
