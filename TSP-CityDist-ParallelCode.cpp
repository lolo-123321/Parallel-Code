#include <iostream> 
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <numeric> 

using namespace std;

const int NUM_CITIES = 12; // Increase this for more computation, e.g., 20-22 cities
vector<pair<int, int>> cities(NUM_CITIES);

// Function to calculate the distance between two cities
double distance(pair<int, int> city1, pair<int, int> city2)
{
    return sqrt(pow(city1.first - city2.first, 2) + pow(city1.second - city2.second, 2));
}

// Function to calculate the total distance of a tour (a specific route) with nested loops
double calculateRouteDistance(const vector<int> &route, const vector<pair<int, int>> &cities)
{
    double total_distance = 0.0;

    // Outer loop for visiting every city in the route
    for (int i = 0; i < route.size() - 1; ++i)
    {
        // Inner loop to simulate repeated recalculations (nested complexity)
        for (int j = 0; j < route.size(); ++j)
        {
            // Calculating the distance between city 'i' and city 'j' repeatedly
            total_distance += distance(cities[route[i]], cities[route[j]]);
        }
    }

    // Add the distance to return to the starting city (repeated with another loop)
    for (int k = 0; k < route.size(); ++k)
    {
        total_distance += distance(cities[route.back()], cities[route[k]]);
    }

    return total_distance;
}

// Brute-force function to solve TSP with nested loops
double solveTSP(const vector<pair<int, int>> &cities)
{
    vector<int> city_indices(cities.size());
    iota(city_indices.begin(), city_indices.end(), 0); 

    double min_distance = numeric_limits<double>::max();

    // Generate all permutations of the cities and calculate the total distance
    do
    {
        // Start the timer for this section
        auto section_start = chrono::high_resolution_clock::now();

        // Nested loop to add complexity while calculating the route distance
        for (int m = 0; m < 1000; ++m)
        { // Outer loop for repeated distance calculations
            double current_distance = calculateRouteDistance(city_indices, cities);
            if (current_distance < min_distance)
            {
                min_distance = current_distance;
            }
        }

        // End the timer for this section
        auto section_end = chrono::high_resolution_clock::now();
        chrono::duration<double> section_diff = section_end - section_start;
        cout << "Section execution time (1000 distance calculations): " << section_diff.count() << " seconds" << endl;

    } while (next_permutation(city_indices.begin() + 1, city_indices.end())); // Skip the first city (start point)

    return min_distance;
}

int main()
{
    //  city coordinates (randomly generated)
    srand(time(0));
    for (int i = 0; i < NUM_CITIES; ++i)
    {
        cities[i] = {rand() % 1000, rand() % 1000}; // Cities are in a 1000x1000 grid
    }

    // Print cities
    cout << "Cities (X, Y): " << endl;
    for (const auto &city : cities)
    {
        cout << "(" << city.first << ", " << city.second << ")" << endl;
    }

    auto start = chrono::high_resolution_clock::now();

    // Solve TSP
    double min_distance = solveTSP(cities);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;

    cout << "Minimum Distance: " << min_distance << endl;
    cout << "Total Execution Time: " << diff.count() << " seconds" << endl;

    return 0;
}
