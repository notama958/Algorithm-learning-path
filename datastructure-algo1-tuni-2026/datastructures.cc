// Datastructures.cc

#include "datastructures.hh"

#include <cmath>
#include <random>
#include <queue>
#include <unordered_set>
#include <iostream>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

Datastructures::Datastructures()
{
    // Write any initialization you need here

}

Datastructures::~Datastructures()
{
    // Write any cleanup you need here

}

/**
 * Adds a beacon to the data structure.
 * 
 * @param id BeaconID
 * @param name Name of the beacon
 * @param xy Coordinate location
 * @param color Color of the beacon 
 * @return true if beacon was added successfully, false if beacon with given ID already exists
 * 
 * Time complexity: O(1) average case (hash map lookup + vector push_back)
 * Space complexity: O(1) amortized
 */
bool Datastructures::add_beacon(BeaconID id, const Name& name, Coord xy, Color color)
{
    // Check for BeaconID existence - O(1) average
    if (m_beaconIdIndices.find(id) != m_beaconIdIndices.end()) {
        return false;
    }
    // Add new beacon to the pool
    m_beacons.emplace_back(std::move(Beacon(id, name, xy, color)));
    // Add save beaconId and its index in the lookup table
    m_beaconIdIndices.emplace(id, m_beacons.size() - 1);
    // Cache the min and max beacon brightness
    updateBrightnessCache(id, color);
    // Reset sorting
    m_prevSortType = SortType::NONE;
    return true;
}

/**
 * @brief Returns the total number of beacons in the system
 * 
 * @return int Current count of all beacons
 * 
 * @complexity O(1) (vector::size())
 */
int Datastructures::beacon_count()
{
    return m_beacons.size();
}

/**
 * @brief Removes all beacons and resets the data structure to empty state
 * 
 * Clears all beacon data, lightbeam connections, and cached information.
 * After calling this function, the system behaves as if newly created.
 * 
 * @post beacon_count() == 0
 * @post all_beacons() returns empty vector
 * @post All ID lookups return sentinel values (NO_NAME, NO_COORD, NO_COLOR)
 * 
 * @complexity O(n) where n is the number of beacons (clearing all containers)
 */
void Datastructures::clear_beacons()
{
    // ===== PRIMARY STORAGE =====
    m_beacons.clear();
    m_beaconIdIndices.clear();
    
    // ===== SORTING CACHE SYSTEM =====
    m_sortedIndices.clear();
    m_prevSortType = SortType::NONE;
    
    // ===== BRIGHTNESS CACHE =====
    m_minBeaconBrightness.reset();
    m_maxBeaconBrightness.reset();
    
    // ===== COLOR CALCULATION CACHE =====
    m_colorCache.clear();
}

/**
 * @brief Retrieves all beacon IDs currently in the system
 * 
 * Returns a list of all beacon ID. The order is implementation-defined
 * and may vary between calls unless the data structure hasn't been modified.
 * 
 * @return std::vector<BeaconID> Vector containing all beacon IDs
 * 
 * @note This operation is not included in performance tests
 * @complexity O(n) where n is beacon_count()
 * @exception Strong exception safety (if allocation fails, structure unchanged)
 */
std::vector<BeaconID> Datastructures::all_beacons()
{
    std::vector<BeaconID> beacons;
    for(const auto& beacon : m_beacons) {
        beacons.emplace_back(beacon.id);
    }
    return beacons;
}

/**
 * @brief Retrieves the name of a beacon by its ID
 * 
 * @param id BeaconID to look up
 * 
 * @return Name of the beacon if found, NO_NAME otherwise
 * 
 * @note This operation is called frequently by the main program
 * @complexity O(1) average case (hash map lookup)
 */
Name Datastructures::get_name(BeaconID id)
{
    const auto foundBeacon = m_beaconIdIndices.find(id);
    if(foundBeacon == m_beaconIdIndices.end()) {
        return NO_NAME;
    }
    return m_beacons[foundBeacon->second].name;
}

/**
 * @brief Retrieves the coordinates of a beacon by its ID
 * 
 * @param id BeaconID to look up
 * 
 * @return Coord of the beacon if found, NO_COORD otherwise
 * 
 * @note This operation is called frequently by the main program
 * @complexity O(1) average case (hash map lookup)
 */
Coord Datastructures::get_coordinates(BeaconID id)
{
    const auto foundBeacon = m_beaconIdIndices.find(id);
    if(foundBeacon == m_beaconIdIndices.end()) {
        return NO_COORD;
    }
    return m_beacons[foundBeacon->second].coord;
}

/**
 * @brief Retrieves the color of a beacon by its ID
 * 
 * Returns the beacon's intrinsic color, not considering any lightbeam effects.
 * 
 * @param id Beacon identifier to look up
 * 
 * @return Color of the beacon if found, NO_COLOR otherwise
 * 
 * @note This operation is called frequently by the main program
 * @complexity O(1) average case (hash map lookup)
 */
Color Datastructures::get_color(BeaconID id)
{
    const auto foundBeacon = m_beaconIdIndices.find(id);
    if(foundBeacon == m_beaconIdIndices.end()) {
        return NO_COLOR;
    }
    return m_beacons[foundBeacon->second].color;
}

/**
 * @brief Returns beacon IDs sorted by name in alphabetical order
 * 
 * Uses cached sorted indices to avoid repeated sorting when the
 * same sort order is requested multiple times.
 * 
 * @return std::vector<BeaconID> Beacon IDs in alphabetical order by name
 * 
 * @complexity O(n log n) for first call, O(n) for subsequent calls with same sort
 */
std::vector<BeaconID> Datastructures::beacons_alphabetically()
{
    if(SortType::BY_NAME_ALPHA != m_prevSortType) {
        sortBeacons(m_sortedIndices, [](const Beacon& a, const Beacon& b){
            return a.name < b.name;
        });
        m_prevSortType = SortType::BY_NAME_ALPHA;

    }    
    std::vector<BeaconID> result;
    result.reserve(m_sortedIndices.size());
    
    for (size_t index : m_sortedIndices) {
        result.push_back(m_beacons[index].id);
    }
    
    return result;
}

/**
 * @brief Returns beacon IDs sorted by increasing brightness
 * 
 * Brightness is calculated from color using luminance formula.
 * Caches the sorted order for efficiency.
 * 
 * @return std::vector<BeaconID> Beacon IDs ordered by increasing brightness
 * 
 * @complexity O(n log n) for first call, O(n) for subsequent calls
 */
std::vector<BeaconID> Datastructures::beacons_brightness_increasing()
{
    if(SortType::BY_BRIGHTNESS_INC != m_prevSortType ) {
        sortBeacons(m_sortedIndices, [](const Beacon& a, const Beacon& b){
            return calculateBrightness(a.color) < calculateBrightness(b.color);
        });
        m_prevSortType = SortType::BY_BRIGHTNESS_INC;
    }    
    std::vector<BeaconID> result;
    result.reserve(m_sortedIndices.size());
    
    for (size_t index : m_sortedIndices) {
        result.push_back(m_beacons[index].id);
    }
    
    return result;
}

/**
 * @brief Returns the beacon ID with minimum brightness
 * 
 * Maintains cached minimum brightness value for O(1) retrieval.
 * Cache is updated when beacons are added or colors changed.
 * 
 * @return BeaconID of the beacon with minimum brightness, or NO_BEACON if none
 * 
 * @complexity O(1) average case (cached value)
 */
BeaconID Datastructures::min_brightness()
{
    if (!m_minBeaconBrightness.has_value()) {
        if (m_beacons.empty()) {
            return NO_BEACON;
        }
        
        auto min_value = std::numeric_limits<int>::max();
        BeaconID min_id = NO_BEACON;
        
        for (const auto& beacon : m_beacons) {
            float brightness = calculateBrightness(beacon.color);
            if (brightness < min_value) {
                min_value = brightness;
                min_id = beacon.id;
            }
        }
        
        m_minBeaconBrightness = {min_id, min_value};
    }
    
    return m_minBeaconBrightness->first;
}

/**
 * @brief Returns the beacon ID with maximum brightness
 * 
 * Maintains cached maximum brightness value for O(1) retrieval.
 * Cache is updated when beacons are added or colors changed.
 * 
 * @return BeaconID of the beacon with maximum brightness, or NO_BEACON if none
 * 
 * @complexity O(1) average case (cached value)
 */
BeaconID Datastructures::max_brightness()
{
    if (!m_maxBeaconBrightness.has_value()) {
        if (m_beacons.empty()) {
            return NO_BEACON;
        }
        
        int max_value = 0;
        BeaconID max_id = NO_BEACON;
        
        for (const auto& beacon : m_beacons) {
            float brightness = calculateBrightness(beacon.color);
            if (brightness > max_value) {
                max_value = brightness;
                max_id = beacon.id;
            }
        }
        
        m_maxBeaconBrightness = {max_id, max_value};
    }
    
    return m_maxBeaconBrightness->first;
}

/**
 * @brief Finds all beacons with the specified name
 * 
 * Optimizes search using binary search if the data is already sorted by name.
 * Otherwise performs linear search.
 * 
 * @param name Name to search for
 * @return std::vector<BeaconID> All beacon IDs with matching name
 * 
 * @complexity O(log n + k) if sorted by name, O(n) otherwise
 *            where n = total beacons, k = matching beacons
 */
std::vector<BeaconID> Datastructures::find_beacons(Name const& name)
{
    std::vector<BeaconID> result;
    // If we're already sorted by name, use binary search
    if (m_prevSortType == SortType::BY_NAME_ALPHA) {
        const auto& sortedIndices = m_sortedIndices;
        // Find first occurrence
        auto it = std::lower_bound(sortedIndices.begin(), sortedIndices.end(), name,
            [this](size_t index, const std::string& targetName) {
            return m_beacons[index].name < targetName;
            });
        
        // Collect all matches starting from found position
        while (it != sortedIndices.end() && 
               m_beacons[*it].name == name) {
            result.emplace_back(m_beacons[*it].id);
            ++it;
        }

    } else {
        // Linear search O(n)
        for(const auto& beacon : m_beacons) {
            if (beacon.name == name) {
                result.push_back(beacon.id);
            }
        }
        std::sort(result.begin(), result.end());
    }
    return result;
}

/**
 * @brief Changes the name of an existing beacon
 * 
 * Updates the beacon's name and invalidates any cached sorts
 * that depend on name ordering.
 * 
 * @param id Beacon ID to modify
 * @param newname New name for the beacon
 * @return true if beacon was found and name changed
 * @return false if no beacon with given ID exists
 * 
 * @post If successful, get_name(id) returns newname
 * @post Any cached name-based sorts are invalidated
 * 
 * @complexity O(1) average case
 */
bool Datastructures::change_beacon_name(BeaconID id, const Name& newname)
{
    auto foundBeacon = m_beaconIdIndices.find(id);
    if(foundBeacon == m_beaconIdIndices.end()) {
        return false;
    }
    m_beacons[foundBeacon->second].name = newname;
    // Invalidate cached sorts that depend on name
    if (m_prevSortType == SortType::BY_NAME_ALPHA) {
        m_prevSortType = SortType::NONE;
    }
    return true;
}


/**
 * @brief Adds a lightbeam from source beacon to target beacon
 * 
 * @param sourceid Source beacon ID
 * @param targetid Target beacon ID
 * @return true if lightbeam added successfully
 * @return false if beacons don't exist, source already has beam, or creates cycle
 * 
 * @complexity O(n) worst case (cycle detection), O(1) average
 */
bool Datastructures::add_lightbeam(BeaconID sourceid, BeaconID targetid)
{
    // Check if both beacons exist
    auto source_it = m_beaconIdIndices.find(sourceid);
    auto target_it = m_beaconIdIndices.find(targetid);
    
    if (source_it == m_beaconIdIndices.end() || 
        target_it == m_beaconIdIndices.end()) {
        return false;
    }

    Beacon& source = m_beacons[source_it->second];
    Beacon& target = m_beacons[target_it->second];
    
    // Check if source already has a target
    if (source.targetId != NO_BEACON) {
        return false;
    }

    // Check final targetId of the target Beacon
    // might be in separate function
    auto targetBeacon = m_beaconIdIndices.find(targetid);
    if(targetBeacon == m_beaconIdIndices.end()) {
        return false;
    }

    // Check if it is cyclic before adding lightbeam 
    BeaconID current = targetid;
    while (current != NO_BEACON) {
        if (current == sourceid) {
            return false;  // Cycle detected
        }
        current = m_beacons[m_beaconIdIndices[current]].targetId;
    }

    // Add the lightbeam
    source.targetId = targetid;
    target.sourceIds.push_back(sourceid);
    // Reset calculated at m_colorCache
    m_colorCache.erase(targetid);
    return true;
}

/**
 * @brief Gets all direct light sources pointing to a beacon
 * 
 * @param id Beacon ID to query
 * @return std::vector<BeaconID> List of beacon IDs that point to this beacon
 * 
 * @note Returns empty vector if no light sources or beacon doesn't exist
 * @complexity O(slogs) where s is number of sources (sorting descending)
 */
std::vector<BeaconID> Datastructures::get_lightsources(BeaconID id)
{
    if(NO_BEACON == id) {
        return { NO_BEACON };
    }
    auto it = m_beaconIdIndices.find(id);
    if (it == m_beaconIdIndices.end()) {
        return { NO_BEACON };  // Beacon not found - NO_BEACON
    }
    
    // Return a copy of the source IDs (already stored in Beacon)
    std::vector<BeaconID> result = m_beacons[it->second].sourceIds;
    
    // Sort in descending order as required
    std::sort(result.begin(), result.end(), std::less<BeaconID>());
    
    return result;
}

/**
 * @brief Gets the path of lightbeams starting from a beacon
 * 
 * Follows the chain of outbeams until reaching a beacon with no outbeam.
 * 
 * @param id Starting beacon ID
 * @return std::vector<BeaconID> Path from start to end, or {NO_BEACON} if invalid
 * 
 * @complexity O(k) where k is path length
 */
std::vector<BeaconID> Datastructures::path_outbeam(BeaconID id)
{
    auto it = m_beaconIdIndices.find(id);
    if (it == m_beaconIdIndices.end()) {
        return { NO_BEACON };
    }

    std::vector<BeaconID> path;
    BeaconID current = id;
    
    // Follow the chain
    while (current != NO_BEACON) {
        path.push_back(current);
        
        auto beacon_it = m_beaconIdIndices.find(current);
        if (beacon_it == m_beaconIdIndices.end()) {
            break;  // Should not happen if data is consistent
        }
        
        current = m_beacons[beacon_it->second].targetId;
    }
    
    return path;    
}

/**
 * @brief Finds the longest incoming beam path to a beacon
 * 
 * Performs DFS on the incoming beam graph to find the longest path
 * ending at the given beacon.
 * 
 * @param id Target beacon ID
 * @return std::vector<BeaconID> Longest path from a source to this beacon
 * 
 * @complexity O(n + e) where n = nodes, e = edges in the subgraph
 */
std::vector<BeaconID> Datastructures::path_inbeam_longest(BeaconID id)
{
    auto sourceBeacon = m_beaconIdIndices.find(id);
    if(sourceBeacon == m_beaconIdIndices.end()) {
        return { NO_BEACON };
    }
    std::vector<BeaconID> maxPath;
    std::vector<BeaconID> tracePath;
    dfs_find_longest(tracePath, m_beacons[sourceBeacon->second], maxPath);
    std::reverse(maxPath.begin(), maxPath.end());
    return maxPath;
}

/**
 * @brief Recursively checking every paths and keep track to possible maximum path found
 * 
 * @param path Current exploring path
 * @param node starting beacon node
 * @return std::vector<BeaconID> reference to the longest path from a source to this beacon
 * 
 * @complexity O(n + e) where n = nodes, e = edges in the subgraph
 */
void Datastructures::dfs_find_longest(std::vector<BeaconID>& path, const Beacon& node, std::vector<BeaconID>& maxPath) 
{
    path.push_back(node.id);
    // Base case: leaf node (no incoming beams)
    if(node.sourceIds.empty()) {
        if(path.size() > maxPath.size()) {
            maxPath = path; // Update longest path found
        }
        return; // Returns to previous level
    } 
    const auto& sourceIds = node.sourceIds;
    for(const auto& sourceId : sourceIds) {
        const auto& sourceNode = m_beacons[m_beaconIdIndices[sourceId]];
        dfs_find_longest(path, sourceNode, maxPath); // recursively check every branch till we reach root node of the branch
        if(!path.empty()) path.pop_back(); // recusively clean up branch to explore another node (branch)
    }
}

/**
 * @brief Calculates total color of a beacon including all incoming beams
 * 
 * Color is calculated as average of beacon's own color and all
 * incoming beam colors (recursively calculated).
 * 
 * @param id Beacon ID
 * @return Color Total calculated color, or NO_COLOR if beacon doesn't exist
 * 
 * @complexity O(n + e) where n = nodes, e = edges in the subgraph
 */
Color Datastructures::total_color(BeaconID id)
{
    auto sourceBeacon = m_beaconIdIndices.find(id);
    if(sourceBeacon == m_beaconIdIndices.end()) {
        return NO_COLOR;
    }
    // Avoid recalculate color if there is no changes
    if (m_colorCache.find(id) != m_colorCache.end()) {
        return m_colorCache[id];
    }

    // Stack to collect colors during DFS
    std::vector<Color> colorStack;

    // Start recursive calculation
    Beacon& start_beacon = m_beacons[sourceBeacon->second];
    auto totalColor = findColor(colorStack, start_beacon);
    // Save the calculated color to cache
    m_colorCache[id] = totalColor; 
    return totalColor;
}

/**
 * @brief Recursive helper that uses stack to collect colors in LIFO order
 * 
 * Algorithm:
 * 1. For each source beacon, recursively calculate its color
 * 2. Push each calculated color onto the stack
 * 3. After processing all sources, pop colors from stack (most recent first)
 * 4. Average popped colors with current beacon's color
 * 
 * This ensures we process the DEEPEST/MOST RECENT colors first (LIFO).
 * 
 * @param color_stack Stack to store intermediate color results
 * @param beacon Current beacon being processed
 * @return Color Calculated average color for this beacon
 */
Color Datastructures::findColor(std::vector<Color>& colorStack, const Beacon& beacon)
{
    const auto& sourceIds = beacon.sourceIds;
    if(sourceIds.empty()) {
        return beacon.color;
    }
    // Process each source
    for(const auto& sourceId : sourceIds) {
        const auto& beacon = m_beacons[m_beaconIdIndices[sourceId]];
        auto beaconColor = findColor(colorStack, beacon);
        colorStack.push_back(beaconColor);
    }
    // Calculate the total color of this beacon
    // by pop out number of sourceIds from colorStack and calculate the avg of all sourceIds
    Color total;
    total.r= 0;
    total.g= 0;
    total.b= 0;
    auto num = sourceIds.size();
    for (size_t i = 0; i < num && !colorStack.empty(); ++i) {
        const auto& tbeacon = colorStack.back();
        total.r += tbeacon.r;
        total.g += tbeacon.g;
        total.b += tbeacon.b;
        colorStack.pop_back();
    }
    // Calculate the avg between target beacon color and total
    total.r = std::round((total.r + beacon.color.r) / ( num + 1 ));  
    total.g = std::round((total.g + beacon.color.g) / ( num + 1 ));  
    total.b = std::round((total.b + beacon.color.b) / ( num + 1 ));  
    return total; 
}

/**
 * @brief Adds a fibre connection between two xpoints with given cost
 * 
 * @param xpoint1 First endpoint coordinate
 * @param xpoint2 Second endpoint coordinate  
 * @param cost Connection cost (distance/weight)
 * @return true if fibre was added, false if already exists or invalid
 * 
 * @complexity O(1) average, O(n) worst-case (hash collision)
 * @rationale Two hash map lookups and insertions, all O(1) average
 */
bool Datastructures::add_fibre(Coord xpoint1, Coord xpoint2, Cost cost)
{
    // Check for self-loop (point connected to itself)
    if(xpoint1 == xpoint2) {
        return false;
    }
    // Check if edge already exists in either direction
    // Using count() is O(1) average for unordered_map
    // Check if xpoint1 exists
    auto it1 = m_fibreGraph.find(xpoint1);
    if (it1 != m_fibreGraph.end() && it1->second.count(xpoint2) > 0) {
        return false;
    }

    // Check if xpoint2 exists
    auto it2 = m_fibreGraph.find(xpoint2);
    if (it2 != m_fibreGraph.end() && it2->second.count(xpoint1) > 0) {
        return false;
    }
    // Add edge in both directions (undirected graph)
    // operator[] creates entry if not exists: O(1) average
    m_fibreGraph[xpoint1][xpoint2] = cost;
    m_fibreGraph[xpoint2][xpoint1] = cost;
    return true;
}

/**
 * @brief Returns all xpoint coordinates in the fibre network
 * 
 * @return vector<Coord> List of all xpoint coordinates
 * 
 * @complexity O(n) where n = number of vertices (xpoints)
 * @rationale Must iterate through all keys in the outer map
 */
std::vector<Coord> Datastructures::all_xpoints()
{
    std::vector<Coord> xpoints;
    // Pre-allocate memory for efficiency
    xpoints.reserve(m_fibreGraph.size());
    // Iterate through all vertices: O(n)
    for (const auto& [coord, neighbors] : m_fibreGraph) {
        xpoints.push_back(coord);
    }
    return xpoints;
}

/**
 * @brief Returns all fibres originating from a specific xpoint
 * 
 * @param xpoint Source coordinate to query
 * @return vector<pair<Coord, Cost>> List of (destination, cost) pairs
 * 
 * @complexity O(k) where k = degree of xpoint (number of neighbors)
 * @rationale Single hash lookup O(1) + iteration over neighbor map O(k)
 */
std::vector<std::pair<Coord, Cost>> Datastructures::get_fibres_from(Coord xpoint)
{   
    std::vector<std::pair<Coord, Cost>> result;
    // Find vertex in graph: O(1) average
    const auto& coordIt = m_fibreGraph.find(xpoint);

    if(coordIt != m_fibreGraph.end()) {
        const auto& coords = coordIt->second;
        // Reserve space for all neighbors
        result.reserve(coords.size());
        // Copy all edges: O(k) where k = degree
        for(auto coord : coords) {
            result.push_back(coord);
        }
    }
    return result;
}

/**
 * @brief Returns all unique fibre connections in the network
 * 
 * @return vector<pair<Coord, Coord>> List of unique (coord1, coord2) pairs
 * 
 * @complexity O(n + e) where n = vertices, e = edges
 * @rationale Iterates through all vertices and their edges
 */
std::vector<std::pair<Coord, Coord> > Datastructures::all_fibres()
{
    std::vector<std::pair<Coord, Coord>> fibres;

    for (auto& [coord1, neighbors] : m_fibreGraph) {
        for (auto& [coord2, cost] : neighbors) {
            // Compare coordinates to include each undirected edge only once
            if (coord1 < coord2) {
                fibres.push_back({coord1, coord2});
            }
        }
    }
    return fibres;
}


/**
 * @brief Removes a fibre connection between two xpoints
 * 
 * @param xpoint1 First endpoint
 * @param xpoint2 Second endpoint
 * @return true if fibre existed and was removed, false otherwise
 * 
 * @complexity O(1) average for lookups, O(k) worst-case for hash collision
 * @rationale Constant time hash lookups and map deletions
 */
bool Datastructures::remove_fibre(Coord xpoint1, Coord xpoint2)
{
    // Find first vertex and check if edge exists: O(1) average
    auto xpoint1Edges = m_fibreGraph.find(xpoint1);
    bool edgeExisted = xpoint1Edges != m_fibreGraph.end() && 
        xpoint1Edges->second.find(xpoint2) != xpoint1Edges->second.end();
    if(!edgeExisted) {
        return false;
    }
    // Remove edge from both adjacency lists: O(1) average
    xpoint1Edges->second.erase(xpoint2);
    auto xpoint2Edges = m_fibreGraph.find(xpoint2); 
    xpoint2Edges->second.erase(xpoint1);

    // Clean up vertices with no remaining edges (optional optimization)
    if (xpoint1Edges->second.empty()) {
        m_fibreGraph.erase(xpoint1);
    }
    if (xpoint2Edges->second.empty()) {
        m_fibreGraph.erase(xpoint2);
    }
    return true;
}

/**
 * @brief Removes all fibres and xpoints from the network
 * 
 * @complexity O(n) where n = total size of all data structures
 * @rationale Clearing unordered_map destructs all elements
 */
void Datastructures::clear_fibres()
{
    m_fibreGraph.clear();
}

/**
 * @brief Finds any route between two xpoints using BFS
 * 
 * Uses breadth-first search to find a path between fromxpoint and toxpoint.
 * Returns the first found path with accumulated costs at each step.
 * If no path exists or points are invalid, returns empty vector.
 * 
 * @param fromxpoint Starting coordinate
 * @param toxpoint Target coordinate
 * @return vector of (coordinate, accumulated_cost) pairs along the path
 */
std::vector<std::pair<Coord, Cost> > Datastructures::route_any(Coord fromxpoint, Coord toxpoint)
{
    // Check if both endpoints exist in the graph: O(1) average
    auto from = m_fibreGraph.find(fromxpoint);
    auto to = m_fibreGraph.find(toxpoint);
    auto end = m_fibreGraph.end();
    if(from == end || to == end) {
        return {}; // One or both points don't exist
    }

    std::queue<const Coord*> q; // Queue for BFS frontier
    std::unordered_map<std::size_t, std::pair<const Coord*, Cost>> parent; // child -> parent using CoordHash as key and pointer to parent and its cost
    std::unordered_set<std::size_t> visited; // save visited xpoints's hashed coordinate

    // Initialize BFS from starting point
    q.push(&from->first);
    const auto fromxpointHash = m_coordHasher(fromxpoint);
    parent[fromxpointHash] = {nullptr, 0}; // Start has no parent, cost 0
    visited.insert(fromxpointHash); // Mark start as visited

    // BFS main loop: explores nodes level by level
    while(!q.empty()) {

        const Coord* current = q.front();
        q.pop();

        // Check if we reached the target
        if(*current == toxpoint) {
             // Reconstruct path from target back to start
            std::vector<std::pair<Coord, Cost>> result;
            const Coord* node = current;
            // Trace back through parent pointers
            while(node != nullptr) {
                auto parent_it = parent.find(m_coordHasher(*node)); // Find parent of the node
                result.push_back({*node, parent_it->second.second}); // Save the node + its total_cost
                node = parent_it->second.first; // Move to parent
            }
            // Reverse to get start→end order
            std::reverse(result.begin(), result.end());
            return result;
        }
        // Get current node's neighbors from adjacency list 
        auto current_it = m_fibreGraph.find(*current);
        if(current_it == m_fibreGraph.end()) continue; // Should not happen if graph is consistent

        // Explore all neighbors
        for(auto& [neighbor, cost] : current_it->second) {
            auto hashedCoord = m_coordHasher(neighbor);
            // Skip already visited neighbors
            if(visited.find(hashedCoord) == visited.end()) {
                    
                visited.insert(hashedCoord);
                // Calculate accumulated cost to reach this neighbor
                auto current_parent_it = parent.find(m_coordHasher(*current));
                Cost total_cost = current_parent_it->second.second + cost;
                // Record parent and cost for path reconstruction
                parent[hashedCoord] = {current, total_cost};
                // Add neighbor to BFS frontier
                q.push(&neighbor);
            }
        }
    }
    // No path found
    return {};
}

/**
 * @brief Finds route with fewest intermediate xpoints (shortest hop count)
 * 
 * Uses breadth-first search to find the path with minimum number of
 * intermediate points between fromxpoint and toxpoint.
 * BFS guarantees finding the shortest path in terms of hop count.
 * Returns empty vector if no path exists or points are invalid.
 * 
 * @param fromxpoint Starting coordinate
 * @param toxpoint Target coordinate  
 * @return vector of (coordinate, accumulated_cost) pairs along shortest path
 */
std::vector<std::pair<Coord, Cost>> Datastructures::route_least_xpoints(Coord fromxpoint, Coord toxpoint)
{
    // Check if both endpoints exist in the graph: O(1) average
    auto from = m_fibreGraph.find(fromxpoint);
    auto to = m_fibreGraph.find(toxpoint);
    auto end = m_fibreGraph.end();
    if(from == end || to == end) {
        return {};  // One or both points don't exist
    }

    std::queue<const Coord*> q; // Queue for BFS frontier
    std::unordered_map<std::size_t, std::pair<const Coord*, Cost>> parent; // child -> parent using CoordHash as key and pointer to parent and its cost
    std::unordered_set<std::size_t> visited; // save visited xpoints's hashed coordinate

    // Initialize BFS from starting point
    q.push(&from->first);
    const auto fromxpointHash = m_coordHasher(fromxpoint);
    parent[fromxpointHash] = {nullptr, 0}; // Start has no parent, cost 0
    visited.insert(fromxpointHash);  // Mark start as visited

    // BFS main loop: explores nodes level by level
    while(!q.empty()) {
        const Coord* current = q.front();
        q.pop();

        // Get current node from graph (should exist)
        auto current_it = m_fibreGraph.find(*current);
        if(current_it == m_fibreGraph.end()) continue;

        // Get parent info for cost calculation
        auto current_parent_it = parent.find(m_coordHasher(*current));

        // Explore all neighbors of current node
        for(auto& [neighbor, cost] : current_it->second) {
            auto hashedCoord = m_coordHasher(neighbor);

            // Check if we found the target
            if(neighbor == toxpoint) {
                // Found shortest path to target (BFS guarantees shortest hops)
                std::vector<std::pair<Coord, Cost>> result;

                // Add target with accumulated cost
                result.push_back({toxpoint, current_parent_it->second.second + cost});
                
                // Reconstruct path backwards to start
                const Coord* node = current;
                while(node != nullptr) {
                    auto parent_it = parent.find(m_coordHasher(*node));
                    result.push_back({*node, parent_it->second.second});
                    node = parent_it->second.first; // // Move to parent
                }
                std::reverse(result.begin(), result.end());
                return result;
            }
             // If neighbor not visited, add to BFS frontier
            if(visited.find(hashedCoord) == visited.end()) {
                visited.insert(hashedCoord);

                // Calculate accumulated cost to reach this neighbor
                Cost total_cost = current_parent_it->second.second + cost;

                // Record parent and cost for path reconstruction
                parent[hashedCoord] = {current, total_cost};

                // Add to BFS queue for further exploration
                q.push(&neighbor);
            }
        }
    }

    // No path exists between the points
    return {};
}

/**
 * @brief Finds the fastest (minimum cost) route between two xpoints using A* algorithm
 * 
 * Uses Dijkstra's algorithm to find the path with minimum total cost between
 * two xpoints in a weighted undirected graph. Returns the optimal path along
 * with accumulated costs at each step.
 * 
 * @param fromxpoint Starting coordinate
 * @param toxpoint Target coordinate  
 * @return vector of (coordinate, accumulated_cost) pairs along optimal path,
 *         or empty vector if no path exists or endpoints are invalid
 * 
 * @note Time complexity: O((V + E) log V) with binary heap
 * @note Space complexity: O(V) for priority queue and distance tracking
 */
std::vector<std::pair<Coord, Cost>> Datastructures::route_fastest(Coord fromxpoint, Coord toxpoint)
{
    // Validate input: both endpoints must exist in the graph
    auto from = m_fibreGraph.find(fromxpoint);
    auto to = m_fibreGraph.find(toxpoint);
    auto end = m_fibreGraph.end();
    if(from == end || to == end) {
        return {}; // Invalid endpoints
    }

    // Priority queue for open set: processes nodes with smallest known distance first
    std::priority_queue<GraphNode, std::vector<GraphNode>, CompareCost> q;

    // Parent map: for path reconstruction
    // Stores {parent_coord_pointer, actual_cost_to_reach_node}
    std::unordered_map<std::size_t, std::pair<const Coord*, Cost>> parent;
    
    // Distance tracking: {current_best_distance, processed_flag}
    // current_best_distance = shortest known distance from start to node
    // processed_flag = true if node has been expanded (optimal distance found)
    using EstimatedDistance = double;
    std::unordered_map<std::size_t, std::pair<EstimatedDistance, bool>> visited; 

    // Heuristic caching: store Euclidean distances to target for efficiency
    // Static cache persists across function calls for same target
    static std::unordered_map<std::size_t, double> shortest_distance;
    
    // Initialize start node
    const auto fromxpointHash = m_coordHasher(fromxpoint);
    const double fromxpoint_d = 0;
    shortest_distance[fromxpointHash] = 0; // shortest distance to start is 0

    q.push({&from->first, fromxpoint_d}); // distance 0 from start
    parent[fromxpointHash] = {nullptr, fromxpoint_d}; // Start has no parent
    visited.emplace(fromxpointHash, std::make_pair(0, false));// Not process yet

    // Main loop: process nodes in order of increasing estimated_cost
    while(!q.empty()) {
        // Get node with smallest estimated_cost from priority_queue
        const auto [current, estimated_cost] = q.top(); 
        q.pop();

        // Check if we reached the goal
        if(*current == toxpoint) {
            // Reconstruct path from goal back to start
            std::vector<std::pair<Coord, Cost>> result;
            const Coord* node = current;
            // Trace back through parent pointers
            while(node != nullptr) {
                auto parent_it = parent.find(m_coordHasher(*node));
                result.push_back({*node, parent_it->second.second});  // Store coord and g-score
                node = parent_it->second.first;  // Move to parent
            }
            // Reverse to get start→goal order
            std::reverse(result.begin(), result.end());
            return result;
        }

        // Get current node's data
        auto current_it = m_fibreGraph.find(*current);
        auto currentHash = m_coordHasher(*current);
        auto current_parent_it = parent.find(currentHash);

        // Expand current node: examine all neighbors
        for(auto& [neighbor, cost] : current_it->second) {
            auto hashedCoord = m_coordHasher(neighbor);
            auto neighbor_it = m_fibreGraph.find(neighbor);
            auto visitedNode = visited.find(hashedCoord);

            // Dijkstra relaxation: new_distance = current_distance + edge_cost
            double new_estimated_distance = current_parent_it->second.second + cost;
            
            if(visitedNode == visited.end()) {
                // First time seeing this neighbor
                parent.insert({hashedCoord, {current, new_estimated_distance}}); // Store parent and g 
                visited.emplace(hashedCoord, std::make_pair(new_estimated_distance, false)); // Store coordinate and total cost
                // insert to priority queue with estimated cost
                q.push({&neighbor_it->first, new_estimated_distance}); // Add to priority queue
            } else if(auto& processed = visitedNode->second.second ; !processed) {
                // Neighbor is in open set (visited but not processed)
                // Check if this path is better than previously found path
                auto& [estimated_distance, status] = visitedNode->second;
                if(estimated_distance > new_estimated_distance) {
                    // Found a better path to this neighbor
                    parent[hashedCoord] = {current, new_estimated_distance}; 
                    // Note: Should also update visited estimated_cost for consistency
                    visitedNode->second.first = new_estimated_distance; 
                    q.push({&neighbor_it->first, new_estimated_distance}); // Re-queue
                }
            }
            // If neighbor already processed, skip it
        }
        // Current node is well visted and processed
        visited[currentHash].second = true;
    }
    // No path exists between the points
    return {};
}

/**
 * @brief Finds a route that forms a loop starting from a given point
 * 
 * Uses depth-first search (DFS) to find a path that starts from startxpoint
 * and eventually returns to a point already in the path, forming a loop.
 * Returns the entire route that ends in the discovered loop.
 * 
 * @param startxpoint Starting coordinate for the search
 * @return vector<Coord> The complete route that forms a loop,
 *         or empty vector if no loop can be found
 */
std::vector<Coord> Datastructures::route_fibre_cycle(Coord startxpoint)
{
    // Check if starting point exists in the graph
    auto start = m_fibreGraph.find(startxpoint);
    if (start == m_fibreGraph.end()) return {};
    
    // Data structures for DFS:
    // visited: tracks globally visited nodes to avoid infinite loops
    // path: stores the current DFS exploration path from start
    std::unordered_map<std::size_t, bool> visited;
    std::vector<Coord> path;

    // Start recursive DFS from the starting point
    return dfs_find_cycle(startxpoint, nullptr, visited, path);
}

/**
 * @brief Recursive DFS helper to find a loop in the fibre network
 * 
 * Explores the graph depth-first, maintaining the current path.
 * When a node is found that's already in the current path,
 * a loop has been discovered.
 * 
 * @param current Current node being explored
 * @param parent Parent node in DFS tree (to avoid backtracking in undirected graph)
 * @param visited Global visited set to avoid revisiting nodes
 * @param path Current DFS exploration path from start
 * @return vector<Coord> The loop path if found, empty vector otherwise
 */
std::vector<Coord> Datastructures::dfs_find_cycle(
    Coord current, 
    Coord* parent,
    std::unordered_map<std::size_t, bool>& visited,
    std::vector<Coord>& path)
{
    std::size_t current_hash = m_coordHasher(current);
    
    // Cycle detection
    // Check if current node is already in the current path
    // This indicates we've found a loop!
    for (size_t i = 0; i < path.size(); ++i) {
        if (path[i] == current) {
            // Found a cycle! The loop is from where 'current' appears in path
            // to the end of path, then back to 'current'
            std::vector<Coord> cycle;
            // Extract the loop portion: from where current appears to end
            for (size_t j = i; j < path.size(); ++j) {
                cycle.push_back(path[j]);
            }
            // Add current again to close the loop
            cycle.push_back(current); 
            return cycle;
        }
    }
    
    // DFS exploration
    // Mark current node as visited globally
    visited[current_hash] = true;
    path.push_back(current);
    
    // Explore all neighbors of current node
    for (const auto& [neighbor, cost] : m_fibreGraph[current]) {
        // In undirected graphs, skip the edge back to parent
        // This prevents trivial backtracking loops
        if (parent && neighbor == *parent) continue;
        
        // Recursively explore the neighbor
        auto result = dfs_find_cycle(neighbor, &current, visited, path);
        if (!result.empty()) {
            return result;  // If recursion found a cycle, propagate it back up
        }
    }
    
    // Backtracking
    // No cycle found from this branch
    // Remove current node from path before returning to parent
    path.pop_back();
    return {}; // No cycle found in this branch
}

/**
 * @brief Calculates the perceived brightness of a color using weighted formula
 * 
 * Uses the standard luminance formula for RGB colors, which approximates
 * how the human eye perceives brightness. Green contributes most (6×),
 * then red (3×), then blue (1×) based on photopic spectral sensitivity.
 * 
 * @param color RGB color to calculate brightness for
 * @return int Brightness value in range [0, 3060] for 8-bit colors (0-255)
 * @note Maximum for 8-bit RGB: 3×255 + 6×255 + 255 = 3060
 */
int Datastructures::calculateBrightness(const Color& color)
{
    return 3 * color.r + 6 * color.g + color.b;
}

/**
 * @brief Updates brightness caches when a beacon's color changes
 * 
 * Maintains up-to-date minimum and maximum brightness beacons
 * for efficient querying. Called whenever a beacon's color is
 * set or modified.
 * 
 * @param id Beacon ID whose color was updated
 * @param color New color of the beacon
 * 
 * @note Complexity: O(1) per update
 */
void Datastructures::updateBrightnessCache(BeaconID id, Color color) {
    int brightness = calculateBrightness(color);
    
    if (!m_minBeaconBrightness.has_value() || 
        brightness < m_minBeaconBrightness->second) {
        m_minBeaconBrightness = {id, brightness};
    }
    
    if (!m_maxBeaconBrightness.has_value() || 
        brightness > m_maxBeaconBrightness->second) {
        m_maxBeaconBrightness = {id, brightness};
    }
}