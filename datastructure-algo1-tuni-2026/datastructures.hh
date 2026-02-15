// Datastructures.hh

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <source_location>
#include <map>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <optional>
#include <stack>


// Type for beacon IDs
using BeaconID = std::string;
using Name = std::string;

// Return value for cases where required beacon was not found
BeaconID const NO_BEACON= "--NO_BEACON--";

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();

// Return value for cases where name values were not found
Name const NO_NAME = "-- NO_NAME --";

// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};

// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (c1.y < c2.y) { return true; }
    else if (c2.y < c1.y) { return false; }
    else { return c1.x < c2.x; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Type for color (RGB)
struct Color
{
    int r = NO_VALUE;
    int g = NO_VALUE;
    int b = NO_VALUE;
};

// Equality and non-equality comparisons for Colors
inline bool operator==(Color c1, Color c2) { return c1.r == c2.r && c1.g == c2.g && c1.b == c2.b; }
inline bool operator!=(Color c1, Color c2) { return !(c1==c2); }

// Return value for cases where color was not found
Color const NO_COLOR = {NO_VALUE, NO_VALUE, NO_VALUE};

// Type for light transmission cost (used only in the second assignment)
using Cost = int;

// Return value for cases where cost is unknown
Cost const NO_COST = NO_VALUE;

// This exception class is there just so that the user interface can notify
// about operations which are not (yet) implemented
class NotImplemented : public std::exception
{
public:
    explicit NotImplemented(std::string const msg = "",
                            const std::source_location location = std::source_location::current())
        : msg_{}
    {
        std::string funcname = location.function_name();
        if (auto namestart = funcname.find_last_of(':'); namestart != std::string::npos)
        { funcname.erase(0, namestart+1); }
        if (auto nameend = funcname.find_first_of('('); nameend != std::string::npos)
        { funcname.erase(nameend, std::string::npos); }
        msg_ = (!msg.empty() ? msg+" in " : "")+funcname+"()";
    }
    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }
private:
    std::string msg_;
};

// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // A operations

    // Estimate of performance:
    // O(1) average - map lookup + vector push_back
    bool add_beacon(BeaconID id, Name const& name, Coord xy, Color color);

    // Estimate of performance:
    // O(1) - vector::size()
    int beacon_count();

    // Estimate of performance:
    // O(n) - clearing all containers
    void clear_beacons();

    // Estimate of performance:
    // O(n) - iterating through all beacons
    std::vector<BeaconID> all_beacons();

    // Estimate of performance:
    // O(1) average - unordered_map lookup
    Name get_name(BeaconID id);

    // Estimate of performance:
    // O(1) average - unordered_map lookup
    Coord get_coordinates(BeaconID id);

    // Estimate of performance:
    // O(1) average - unordered_map lookup
    Color get_color(BeaconID id);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // Time complexity: O(n log n) for first call or after data modification,
    //                  O(n) for subsequent calls with same sort type
    // Space complexity: O(n) for storing sorted indices
    std::vector<BeaconID> beacons_alphabetically();

    // Estimate of performance:
    // Time complexity: O(n log n) for first call or after data modification,
    //                  O(n) for subsequent calls with same sort type
    // Space complexity: O(n) for storing sorted indices
    std::vector<BeaconID> beacons_brightness_increasing();

    // Estimate of performance:
    // Time complexity: O(1) average case (cached), O(n) worst case (cache miss)
    // Space complexity: O(1)
    BeaconID min_brightness();

    // Estimate of performance:
    // Time complexity: O(1) average case (cached), O(n) worst case (cache miss)
    // Space complexity: O(1)
    BeaconID max_brightness();

    // Estimate of performance:
    // Time complexity: O(log n + k) if data sorted by name, O(n) otherwise
    //                  where n = total beacons, k = number of matches
    // Space complexity: O(k) for result storage
    // where n = total beacons, k = matching beacons
    std::vector<BeaconID> find_beacons(Name const& name);

    // Estimate of performance:
    // - O(1) hash map lookup to find beacon by ID
    // - O(1) to update the name field
    // - O(1) to invalidate relevant caches if needed
    // - No reallocation or data movement required
    bool change_beacon_name(BeaconID id, Name const& newname);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // - Beacon existence check: O(1) hash map lookup
    // - Cycle detection: O(k) where k is chain length (worst case O(n))
    // - Actual beam addition: O(1) (vector push_back amortized)
    bool add_lightbeam(BeaconID sourceid, BeaconID targetid);

    // Estimate of performance:
    // Time complexity: O(s log s) where s = number of sources
    // Space complexity: O(s) for result vector
    std::vector<BeaconID> get_lightsources(BeaconID id);

    // Estimate of performance:
    // Time complexity: O(k) where k = path length. Follows chain until NO_BEACON or visited all (max n steps)
    // Space complexity: O(k) for result vector
    std::vector<BeaconID> path_outbeam(BeaconID id);

    // B operations

    // Estimate of performance:
    // Time complexity: O(n + e) Using DFS traversal where n = nodes, e = edges in subgraph
    // Space complexity: O(n) for recursion stack and path storage
    std::vector<BeaconID> path_inbeam_longest(BeaconID id);

    // Estimate of performance:
    // Time complexity: O(n + e) with memoization, O(2^n) without
    // Space complexity: O(n) for recursion stack and memoization cache
    // - Without memoization: exponential due to repeated calculations
    // - With memoization: each node's color calculated once
    // - Each edge traversed once during recursive calculation
    Color total_color(BeaconID id);

    // Estimate of performance:
    // Time complexity: O(1) average, O(n) worst-case (hash collisions)
    // Space complexity: O(1) additional space
    // - Hash map lookups and insertions, all O(1) average
    // - Rejects self-loops and duplicate edges early
    // - Stores edge in both directions (undirected graph)
    bool add_fibre(Coord xpoint1, Coord xpoint2, Cost cost);

    // Estimate of performance:
    // Time complexity: O(n) where n = number of vertices
    // Space complexity: O(n) for returned vector
    // - Iterates through all vertices in adjacency list
    // - Returns copy of all coordinate keys
    std::vector<Coord> all_xpoints();

    // Estimate of performance:
    // Time complexity: O(k) where k = vertex degree
    // Space complexity: O(k) for output vector
    // - Constant time vertex lookup
    // - Linear iteration through neighbors
    std::vector<std::pair<Coord, Cost>> get_fibres_from(Coord xpoint);

    // Estimate of performance:
    // Time complexity: O(n + e) where n = vertices, e = edges
    // Space complexity: O(e) for output vector
    // - Visits each vertex and each edge once
    // - Filters duplicates for undirected representation
    std::vector<std::pair<Coord, Coord>> all_fibres();

    // Estimate of performance:
    // Time complexity: O(1) average, O(n) worst-case
    // Space complexity: O(1) additional space
    // - Constant time edge removal from both ends
    // - Cleans up isolated vertices automatically
    bool remove_fibre(Coord xpoint1, Coord xpoint2);

    // Estimate of performance:
    // Time complexity: O(n) where n = total elements
    // Space complexity: O(1) additional space
    // - Linear cleanup of all graph data
    // - Complete memory deallocation
    void clear_fibres();

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // Time complexity: O(V + E) in worst case, where V = vertices, E = edges
    // Space complexity: O(V) for BFS data structures
    // - BFS visits each vertex at most once: O(V)
    // - Each edge examined at most once: O(E)
    // - Uses queue, parent map, and visited set, each O(V) space
    // - Path reconstruction adds O(P) where P = path length
    std::vector<std::pair<Coord, Cost>> route_any(Coord fromxpoint, Coord toxpoint);

    // C operations

    // Estimate of performance:
    // Time complexity: O(V + E) where V = vertices, E = edges
    // Space complexity: O(V) for BFS data structures
    // - BFS visits each vertex at most once: O(V)
    // - Each edge examined at most once: O(E)
    // - Uses queue, parent map, and visited set, each O(V) space
    // - Returns first found shortest path in terms of hop count
    std::vector<std::pair<Coord, Cost>> route_least_xpoints(Coord fromxpoint, Coord toxpoint);

    // Estimate of performance:
    // Time complexity: O((V + E) log V) where V = vertices, E = edges
    // Space complexity: O(V) for priority queue and distance maps
    // - Dijkstra's algorithm finds minimum cost path in weighted graphs
    // - Each vertex processed once, each edge examined once
    // - Priority queue operations: O(log V) per insertion/removal
    std::vector<std::pair<Coord, Cost>> route_fastest(Coord fromxpoint, Coord toxpoint);

    // Estimate of performance:
    // Time complexity: O(V + E) in worst case, where V = vertices, E = edges
    // Space complexity: O(V) for recursion stack and path storage
    // - Depth-first search that explores all possible paths from start
    // - Returns first found loop (if any) in the exploration order
    // - Each node visited at most once due to global visited tracking
    // - Path vector stores current exploration route (max depth = V)
    std::vector<Coord> route_fibre_cycle(Coord startxpoint);

private:
    /**
     * @brief Main data structure representing a beacon in the network
     * 
     * Stores all information about a single beacon including its connections,
     * position, color, and metadata. Designed for efficient traversal and
     * relationship tracking in the beam network.
     */
    struct Beacon {
        BeaconID id;                          ///< Unique identifier for the beacon
        BeaconID targetId{NO_BEACON};         ///< Single outgoing beam target (if any)
        std::vector<BeaconID> sourceIds;      ///< All incoming beams from other beacons
        Name name;                            ///< Human-readable name of the beacon
        Coord coord;                          ///< 2D position coordinates (x, y)
        Color color;                          ///< RGB color of the beacon
        
        Beacon(BeaconID bId, Name bName, Coord bCoord, Color bColor) 
            : id(bId), 
              name(std::move(bName)),   // Move string to avoid copy
              coord(bCoord), 
              color(bColor) 
        {}
    };

    // ===== BEACON STORAGE AND INDEXING =====

    /// Primary storage container for all beacons in insertion order
    /// Uses vector for contiguous memory and cache efficiency
    std::vector<Beacon> m_beacons; 
    /// Fast ID lookup: maps BeaconID to index in m_beacons vector
    /// O(1) average case lookup, prevents duplicate IDs                        
    std::unordered_map<BeaconID, size_t> m_beaconIdIndices; 

    // ===== SORTING CACHE SYSTEM =====

    /**
     * @brief Enumeration of supported sorting criteria
     * 
     * Used to cache sorted indices to avoid repeated sorting operations.
     * The cache is invalidated when beacons are added/removed.
     */
    enum SortType {
        NONE = 0,               ///< No cached sort available
        BY_NAME_ALPHA = 1,      ///< Sorted alphabetically by beacon name
        BY_BRIGHTNESS_INC = 2   ///< Sorted by increasing brightness
    };
    /// Tracks the type of last performed sort for cache validation
    mutable SortType m_prevSortType{SortType::NONE};  
    /// Cached vector of indices into m_beacons sorted by last criteria
    /// Maintains original storage order while providing sorted access
    mutable std::vector<size_t> m_sortedIndices;      
    
    // ===== BRIGHTNESS CACHE =====

    /// Type alias for brightness cache entries: (beacon_id, brightness_value)
    using BeaconBrightness = std::pair<BeaconID, int>;  
    
    /// Cache for beacon with minimum/maximum brightness
    mutable std::optional<BeaconBrightness> m_minBeaconBrightness;
    mutable std::optional<BeaconBrightness> m_maxBeaconBrightness;
    
    // ===== ALGORITHM HELPERS AND CACHES =====

    /// Explores the directed beam graph to find the longest chain of
    /// consecutive incoming beams starting from a given node.
    void dfs_find_longest(std::vector<BeaconID>& path, const Beacon& node, std::vector<BeaconID>& maxPath); 
    /// Calculates the average color of a beacon including all incoming beams
    /// using depth-first traversal. Results are cached to avoid recomputation.
    Color findColor(std::vector<Color>& stack, const Beacon& beacon);
    
    /// Cache for memoized color calculations
    /// Maps BeaconID to its precomputed total_color result
    /// Invalidated when beam relationships change
    std::unordered_map<BeaconID, Color> m_colorCache; 

    // ===== UTILITY FUNCTIONS =====

    /// Updates brightness caches when beacon color changes
    void updateBrightnessCache(BeaconID id, Color color);
    /// Calculates perceived brightness from RGB color
    static int calculateBrightness(const Color& color);

    /**
     * @brief Sorts beacon indices according to comparator
     * 
     * @tparam Compare Comparison function type
     * @param indices Output vector of sorted indices
     * @param comp Comparator function object
     * 
     * @complexity O(n log n) for sorting
     */
    template<typename Compare>
    void sortBeacons(std::vector<size_t>& indices, Compare comp) const {
        indices.resize(m_beacons.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        std::sort(indices.begin(), indices.end(),
            [this, &comp](size_t a, size_t b) { 
                return comp(m_beacons[a], m_beacons[b]);
            });
    }

    // ===== FIBRE NETWORK DATA STRUCTURES =====

    /**
     * @brief Adjacency list representation of undirected fibre graph
     * 
     * Two-level map structure:
     * - Outer map: Coordinate → map of neighbors
     * - Inner map: Neighbor coordinate → edge cost
     * 
     * Provides O(log n) neighbor lookup and O(1) edge cost access.
     * Undirected edges are stored in both directions.
     * It benefits with maintain ordered by Coordinate without doing extra sorting 
     * when get_fibre_from, all_fibres
     */
    std::map<Coord, std::map<Coord, Cost>> m_fibreGraph;

    /// Hash function object for Coord to use in unordered containers
    CoordHash m_coordHasher;
    
    // ===== GRAPH ALGORITHM TYPES =====
    
    /// Node representation for priority queues in graph algorithms
    /// Stores (coordinate pointer, cost/distance) pair
    using GraphNode = std::pair<const Coord*, Cost>;

    /**
     * @brief Comparator for GraphNode in min-heap priority queues
     * 
     * Used by Dijkstra algorithm to process nodes with
     * smallest cost/distance first.
     */
    struct CompareCost {
        bool operator()(const GraphNode& a, const GraphNode& b) {
            return a.second > b.second; 
        }
    };
    
    // ===== CYCLE DETECTION HELPER =====

    /// Explores the fibre network depth-first to find paths that form loops.
    /// Returns the first cycle found starting from the given point.
    std::vector<Coord> dfs_find_cycle(
        Coord current, 
        Coord* parent,
        std::unordered_map<std::size_t, bool>& visited,
        std::vector<Coord>& path);

};

#endif // DATASTRUCTURES_HH
