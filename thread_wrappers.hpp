#include <vector>
#include <functional> 
#include <thread>


template<class F, typename... Args>
void run_function_with_threads(size_t n_items, size_t threads, F f, Args... args) {
    std::vector<std::thread> threads_vector;
    size_t batch_size = (n_items / threads) + 1;
    for (size_t step = 0; step < threads; ++step) {
        size_t start = step * batch_size;
        size_t end = (step + 1) * batch_size;
        if (start >= n_items) {
            continue;
        }
        if (end > n_items) {
            end = n_items;
        }
        threads_vector.push_back(
            std::thread(f,
                        start,
                        end,
                        step,
                        std::ref(args)...
            )
        );
    }
    for (size_t i = 0; i < threads_vector.size(); i++) {
        threads_vector[i].join();
    }
}

template<class F, typename Contraner, typename... Args>
void run_function_with_threads2D(size_t n_items, size_t threads, F f, const Contraner trs_vector, Args... args) {
    std::vector<std::thread> threads_vector;
    size_t batch_size = (n_items / threads) + 1;
    size_t queue = 0;
    size_t prev_i = 0;
    size_t i = 0;
    size_t prev_point = 0;
    size_t point = 0;
    size_t step = 0;
    for (i = 0; i < trs_vector.size(); ++i) {
        for (size_t j = i+1; j < trs_vector.size(); ++j) {
            queue += 1;        
            point += 1;
        }
        if (queue > batch_size) {
            size_t end = i + 1;
            threads_vector.push_back(
                std::thread(f,
                        prev_i,
                        end,
                        step,
                        prev_point,
                        std::ref(args)...
               )
            );
            prev_i = i+1;
            prev_point += queue;
            queue = 0;
            step += 1;
        }
    }
    size_t end = i + 1;
    threads_vector.push_back(
                std::thread(f,
                        prev_i,
                        end,
                        step,
                        prev_point,
                        std::ref(args)...
               )
            );
    for (size_t i = 0; i < threads_vector.size(); i++) {
        threads_vector[i].join();
    }
}

template<class F, typename Contraner, typename... Args>
void run_function_with_results_vector(size_t n_items, size_t threads, F f, Contraner trs_vector, Args... args) {

    std::vector<std::thread> threads_vector;
    size_t batch_size = (n_items / threads) + 1;
    size_t queue = 0;
    size_t prev_i = 0;
    size_t i = 0;
    size_t thread_i = 0;
    for (i = 0; i < trs_vector.size(); ++i) {
        for (size_t j = i+1; j < trs_vector.size(); ++j) {
            queue += 1;        
        }
        if (queue > batch_size) {
            size_t end = i + i;
            threads_vector.push_back(
                std::thread(f,
                            prev_i,
                            end,
                            thread_i,
                            std::ref(args)...
                )
            );
            prev_i = i+1;
            queue = 0;
            thread_i++;
        }
    }

    threads_vector.push_back(
        std::thread(f,
                    prev_i,
                    i,
                    thread_i,
                    std::ref(args)...
        )
    );

    for (size_t i = 0; i < threads_vector.size(); ++i) {
        threads_vector[i].join();
    }
}