#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric> // For iota
#include <omp.h>
#include <vector>

using namespace std;

const int NUM_CITIES =
    12; // Increase this for more computation, e.g., 20-22 cities
// vector<pair<int, int>> cities(NUM_CITIES);

vector<pair<int, int>> cities = {
    {595, 759}, {207, 752}, {3, 925}, {492, 958}, {890, 80}, {562, 85}, {463, 217}, {239, 15}, {716, 698}, {596, 137}, {281, 753}, {323, 719}};

// Function to calculate the eiclidean distance between two cities
double distance(pair<int, int> city1, pair<int, int> city2)
{
    return sqrt(pow(city1.first - city2.first, 2) +
                pow(city1.second - city2.second, 2));
}

// Function to calculate the total distance of a tour (a specific route) with
// nested loops
double calculateRouteDistance(const vector<int> &route,
                              const vector<pair<int, int>> &cities)
{
    double total_distance = 0.0;

#pragma omp parallel for num_threads(12)
    for (int i = 0; i < route.size() - 1; ++i)
    {
        double local_distance = 0.0;

        // Inner loop for additional complexity
        for (int j = 0; j < route.size(); ++j)
        {
            local_distance += distance(cities[route[i]], cities[route[j]]);
        }

#pragma omp atomic
        total_distance += local_distance;
    }

    //  Add the distance to return to the starting city (repeated with another
    //  loop)

    for (int k = 0; k < route.size(); ++k)
    {
        double local_distance = 0.0;

        local_distance += distance(cities[route.back()], cities[route[k]]);

#pragma omp atomic
        total_distance += local_distance;
    }

    return total_distance;
}

// Function to calculate factorial for approximating permutations
long long factorial(int n) { return (n <= 1) ? 1 : n * factorial(n - 1); }

// Brute-force function to solve TSP with nested loops
double solveTSP(const vector<pair<int, int>> &cities)
{
    vector<int> city_indices(cities.size());
    iota(city_indices.begin(), city_indices.end(),
         0); // Fill with 0, 1, 2,... N-1

    double min_distance =
        numeric_limits<double>::max(); // Store the minimum distance found so far
                                       // Generate all permutations of the cities
                                       // and calculate the total distance

#pragma omp parallel for num_threads(12)
    for (int i = 0; i < factorial(city_indices.size() - 1);
         ++i) // Approximation for permutations
    {

        vector<int> local_route = city_indices;

        // Generate the i-th permutation
        size_t idx = i;
        for (size_t j = 1; j < local_route.size(); ++j)
        {
            swap(local_route[j], local_route[j + idx % (local_route.size() - j)]);
            idx /= (local_route.size() - j);
        }

        // Calculate the distance for the current route
        double current_distance = calculateRouteDistance(city_indices, cities);

        // Update the shared min_distance (race condition here)
        if (current_distance < min_distance)
        {
            min_distance = current_distance; // Shared variable accessed unsafely
        }
    }

    return min_distance;
}

int main()
{
    // Predefined city coordinates (randomly generated)

    /*   srand(time(0));
      for (int i = 0; i < NUM_CITIES; ++i)
      {
          cities[i] = {rand() % 1000, rand() % 1000}; // Cities are in a 1000x1000
      grid
      } */
    // omp_set_num_threads(4);

    // Print cities
    cout << "Cities (X, Y): " << endl;
    for (const auto &city : cities)
    {
        cout << "(" << city.first << ", " << city.second << ")" << endl;
    }

    auto start = chrono::high_resolution_clock::now();

    double min_distance = solveTSP(cities);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;

    cout << "Minimum Distance: " << min_distance << endl;
    cout << "Execution Time (solveTSP): " << diff.count() << " seconds" << endl;

    return 0;
}
