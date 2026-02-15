// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)

#include <string>
using std::string;
using std::getline;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::flush;
using std::noskipws;

#include <iomanip>
using std::setfill;
using std::setw;

#include <istream>
using std::istream;

#include <ostream>
using std::ostream;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::istringstream;
using std::ostringstream;
using std::stringstream;

#include <iomanip>
using std::setw;

#include <tuple>
using std::tuple;
using std::make_tuple;
using std::get;

#include <regex>
using std::regex_match;
using std::regex_search;
using std::smatch;
using std::regex;
using std::sregex_token_iterator;

#include <algorithm>
using std::find_if;
using std::find;
using std::reverse;
using std::binary_search;
using std::max_element;
using std::max;
using std::min;
using std::sort;

#include <random>
using std::minstd_rand;
using std::uniform_int_distribution;

#include <vector>
using std::vector;

#include <array>
using std::array;

#include <ctime>
using std::time;

#include <utility>
using std::pair;

#include <cstdlib>
using std::div;

#include <limits>
using std::numeric_limits;

#include <cassert>


#include "mainprogram.hh"

#include "datastructures.hh"

#ifdef GRAPHICAL_GUI
#include "mainwindow.hh"
#endif

string const MainProgram::PROMPT = "> ";

MainProgram::CmdResult MainProgram::cmd_beacon_info(std::ostream& /*output*/, MatchIter begin, MatchIter end)
{
    BeaconID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    return {ResultType::IDLIST, MainProgram::CmdResultIDs{id}};
}

void MainProgram::test_get_functions(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there are no beacons
    {
        BeaconID id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.get_name(id);
        ds_.get_coordinates(id);
        ds_.get_color(id);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_add_beacon(ostream& output, MatchIter begin, MatchIter end)
{
    BeaconID id = *begin++;
    string name = *begin++;
    string xstr = *begin++;
    string ystr = *begin++;
    string rstr = *begin++;
    string gstr = *begin++;
    string bstr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    Coord xy = {convert_string_to<int>(xstr), convert_string_to<int>(ystr)};
    int r = convert_string_to<int>(rstr);
    int g = convert_string_to<int>(gstr);
    int b = convert_string_to<int>(bstr);

    bool success = ds_.add_beacon(id, name, xy, {r, g, b});

    view_dirty = true;
    if (success)
    {
        return {ResultType::IDLIST, MainProgram::CmdResultIDs{id}};
    }
    else
    {
        output << "Adding beacon failed! (false returned)" << std::endl;
        return {ResultType::IDLIST, MainProgram::CmdResultIDs{}};
    }
}

MainProgram::CmdResult MainProgram::cmd_change_name(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    BeaconID id = *begin++;
    string newname = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    bool success = ds_.change_beacon_name(id, newname);

    view_dirty = true;
    if (success)
    {
        return {ResultType::IDLIST, MainProgram::CmdResultIDs{id}};
    }
    else
    {
        output << "Changing beacon name failed! (false returned)" << std::endl;
        return {ResultType::IDLIST, MainProgram::CmdResultIDs{}};
    }
}

void MainProgram::test_change_name(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        auto newname = n_to_name(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.change_beacon_name(id, newname);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_add_lightbeam(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    BeaconID sourceid = *begin++;
    BeaconID targetid = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    bool ok = ds_.add_lightbeam(sourceid, targetid);
    if (ok)
    {
        try
        {
            auto sourcename = ds_.get_name(sourceid);
            auto targetname = ds_.get_name(targetid);
            output << "Added lightbeam: " << sourcename << " -> " << targetname << endl;
        }
        catch (NotImplemented&)
        {
            output << "Added a lightbeam (and get_name() threw NotImplemented)." << endl;
        }
    }
    else
    {
        output << "Adding lightbeam failed!" << endl;
    }

    view_dirty = true;
    return {};
}

MainProgram::CmdResult MainProgram::cmd_add_fibre(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string x1str = *begin++;
    string y1str = *begin++;
    string x2str = *begin++;
    string y2str = *begin++;
    string coststr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    Coord xy1 = {convert_string_to<int>(x1str), convert_string_to<int>(y1str)};
    Coord xy2 = {convert_string_to<int>(x2str), convert_string_to<int>(y2str)};
    Cost cost = convert_string_to<int>(coststr);
    bool ok = ds_.add_fibre(xy1, xy2, cost);
    if (ok)
    {
        output << "Added fibre: (" << xy1.x << "," << xy1.y << ") <-> (" << xy2.x << "," << xy2.y << "), cost " << cost << endl;
    }
    else
    {
        output << "Adding fibre failed!" << endl;
    }

    view_dirty = true;
    return {};
}

MainProgram::CmdResult MainProgram::cmd_remove_fibre(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string x1str = *begin++;
    string y1str = *begin++;
    string x2str = *begin++;
    string y2str = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    Coord xy1 = {convert_string_to<int>(x1str), convert_string_to<int>(y1str)};
    Coord xy2 = {convert_string_to<int>(x2str), convert_string_to<int>(y2str)};
    bool ok = ds_.remove_fibre(xy1, xy2);
    if (ok)
    {
        output << "Removed fibre: (" << xy1.x << "," << xy1.y << ") <-> (" << xy2.x << "," << xy2.y << ")" << endl;
    }
    else
    {
        output << "Removing fibre failed!" << endl;
    }

    view_dirty = true;
    return {};
}

void MainProgram::test_remove_fibre(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id1 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        auto id2 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.remove_fibre(ds_.get_coordinates(id1), ds_.get_coordinates(id2));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_path_outbeam(std::ostream& /*output*/, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    BeaconID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto result = ds_.path_outbeam(id);
    if (result.empty()) { return {ResultType::HIERARCHY, MainProgram::CmdResultIDs{NO_BEACON}}; }
    else { return {ResultType::HIERARCHY, result}; }
}

void MainProgram::test_path_outbeam(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.path_outbeam(id);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_path_inbeam_longest(std::ostream& /*output*/, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    BeaconID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto result = ds_.path_inbeam_longest(id);
    if (result.empty()) { return {ResultType::HIERARCHY, MainProgram::CmdResultIDs{NO_BEACON}}; }
    else { return {ResultType::HIERARCHY, result}; }
}

void MainProgram::test_path_inbeam_longest(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.path_inbeam_longest(id);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_total_color(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    BeaconID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto result = ds_.total_color(id);
    auto name = ds_.get_name(id);
    output << "Total color of " << name << ": (" << result.r << "," << result.g << "," << result.b << ")" << endl;

    return {ResultType::IDLIST, MainProgram::CmdResultIDs{id}};
}

void MainProgram::test_total_color(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.total_color(id);
        watch.stop();
    }
}

void MainProgram::add_random_beacons(unsigned int size, Stopwatch* watchp, Coord min, Coord max)
{
    for (unsigned int i = 0; i < size; ++i)
    {
        string name = n_to_name(random_beacons_added_);
        BeaconID id = n_to_id(random_beacons_added_);

        int x = random<int>(min.x, max.x);
        int y = random<int>(min.y, max.y);
        int r = random<int>(1, 255);
        int g = random<int>(1, 255);
        int b = random<int>(1, 255);

        if (watchp) { watchp->start(); }
        ds_.add_beacon(id, name, {x, y}, {r, g, b});
        if (watchp) { watchp->stop(); }

        // Add random target beacon whose number is smaller, with 80 % probability
        if (random_beacons_added_ > 0 && random(0, 100) < 80)
        {
            BeaconID taxerid = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
            ds_.add_lightbeam(id, taxerid);
        }

        ++random_beacons_added_;
    }
}

MainProgram::CmdResult MainProgram::cmd_random_add(ostream& output, MatchIter begin, MatchIter end)
{
    string sizestr = *begin++;
    string minxstr = *begin++;
    string minystr = *begin++;
    string maxxstr = *begin++;
    string maxystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    unsigned int size = convert_string_to<unsigned int>(sizestr);

    Coord def_min{1, 1};
    Coord def_max{100, 100};
    Coord min = def_min;
    Coord max = def_max;
    if (!minxstr.empty() && !minystr.empty() && !maxxstr.empty() && !maxystr.empty())
    {
        min.x = convert_string_to<unsigned int>(minxstr);
        min.y = convert_string_to<unsigned int>(minystr);
        max.x = convert_string_to<unsigned int>(maxxstr);
        max.y = convert_string_to<unsigned int>(maxystr);
    }
    else
    {
        auto beacons = ds_.all_beacons();
        if (!beacons.empty())
        {
            // Find out bounding box
            min = {numeric_limits<int>::max(), numeric_limits<int>::max()};
            max = {numeric_limits<int>::min(), numeric_limits<int>::min()};
            for (auto const& beacon : beacons)
            {
                auto [x,y] = ds_.get_coordinates(beacon);
                if (x < min.x) { min.x = x; }
                if (y < min.y) { min.y = y; }
                if (x > max.x) { max.x = x; }
                if (y > max.y) { max.y = y; }
            }
        }
    }

    if (min == max)
    {
        min = def_min;
        max = def_max;
    }

    add_random_beacons(size, nullptr, min, max);

    output << "Added: " << size << " beacons." << endl;

    view_dirty = true;

    return {};
}

void MainProgram::test_random_add(Stopwatch& watch)
{
    add_random_beacons(1, &watch);
}

void MainProgram::test_extra_add(Stopwatch& /*watch*/)
{
    // Add a beacon without adding it to the measured time
    add_random_beacons(1, nullptr);
}

MainProgram::CmdResult MainProgram::cmd_beacon_count(ostream& output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    output << "Number of beacons: " << ds_.beacon_count() << endl;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_all_beacons(ostream& output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    auto beacons = ds_.all_beacons();
    if (beacons.empty())
    {
        output << "No beacons!" << endl;
    }

    return {ResultType::IDLIST, beacons};
}

void MainProgram::test_all_beacons(Stopwatch& watch)
{
    watch.start();
    ds_.all_beacons();
    watch.stop();
}

MainProgram::CmdResult MainProgram::cmd_all_xpoints(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    auto xpoints = ds_.all_xpoints();
    if (xpoints.empty())
    {
        output << "No xpoints!" << endl;
    }

    return {ResultType::COORDLIST, xpoints};
}

void MainProgram::test_all_xpoints(Stopwatch& watch)
{
    watch.start();
    ds_.all_xpoints();
    watch.stop();
}

MainProgram::CmdResult MainProgram::cmd_all_fibres(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    auto fibres = ds_.all_fibres();
    if (fibres.empty())
    {
        output << "No fibres!" << endl;
    }
    else
    {
        for (auto& [xy1, xy2] : fibres)
        {
            output << "(" << xy1.x << "," << xy1.y << ") -> (" << xy2.x << "," << xy2.y << ")" << endl;
        }
    }

    return {};
}

MainProgram::CmdResult MainProgram::cmd_lightsources(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto sources = ds_.get_lightsources(id);
    if (sources.empty())
    {
        output << "No lightsources!" << endl;
        return {ResultType::SOURCELIST, CmdResultIDs{id}};
    }

    sources.insert(sources.begin(), id); // Add destination id
    return {ResultType::SOURCELIST, sources};
}

void MainProgram::test_lightsources(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.get_lightsources(id);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_fibres(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string xstr = *begin++;
    string ystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    Coord xy = {convert_string_to<int>(xstr), convert_string_to<int>(ystr)};
    auto fibres = ds_.get_fibres_from(xy);
    if (fibres.empty())
    {
        output << "No fibres!" << endl;
    }

    return {ResultType::FIBRELIST, CmdResultFibres{xy, fibres}};
}

void MainProgram::test_fibres(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.get_fibres_from(ds_.get_coordinates(id));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_random_fibres(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string sizestr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    if (ds_.beacon_count() < 2)
    {
        output << "Not enough beacons to add fibres!" << endl;
        return {};
    }

    unsigned int random_fibres = convert_string_to<unsigned int>(sizestr);

    add_random_fibres(output, random_fibres);

    output << "Added at most " << random_fibres << " random fibres." << endl;

    view_dirty = true;

    return {};
}

void MainProgram::test_random_fibres(Stopwatch& watch)
{
    if (random_beacons_added_ > 0) // Don't do anything if there's no beacons
    {
        for (auto i=0; i<10; ++i)
        {
            auto id1 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
            auto id2 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
            Cost cost = random(0, 100);
            watch.start();
            ds_.add_fibre(ds_.get_coordinates(id1), ds_.get_coordinates(id2), cost);
            watch.stop();
        }
    }
}

MainProgram::CmdResult MainProgram::cmd_random_labyrinth(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string xstr = *begin++;
    string ystr = *begin++;
    string extrastr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int xsize = convert_string_to<int>(xstr);
    int ysize = convert_string_to<int>(ystr);
    int extra = convert_string_to<int>(extrastr);

    create_fibre_labyrinth(output, xsize, ysize, extra);

    output << "Added fibre labyrinth with at most " << extra << " extra routes." << endl;

    view_dirty = true;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_clear_beacons(ostream& output, MatchIter begin, MatchIter end)
{
    assert(begin == end && "Invalid number of parameters");

    ds_.clear_beacons();
    init_primes();

    output << "Cleared all beacons" << endl;

    view_dirty = true;

    return {};
}

std::string MainProgram::print_beacon(BeaconID id, ostream& output)
{
    try
    {
        if (id != NO_BEACON)
        {
            auto name = ds_.get_name(id);
            if (name.empty()) { name = "<empty>"; }
            output << name << ": ";

            auto coord = ds_.get_coordinates(id);
            auto [x, y] = coord;
            output << "pos=";
            if (coord != NO_COORD)
            {
                output << "(";
                if (x != NO_VALUE) { output << x; }
                else {output << "NO_VALUE"; }
                output << ",";
                if (y != NO_VALUE) { output << y; }
                else {output << "NO_VALUE"; }
                output << ")";
            }
            else { output << "NO_COORD"; }
            output << ", ";

            Color color = ds_.get_color(id);
            output <<  "color=";
            if (color != NO_COLOR)
            {
                int brightness = 3*color.r + 6*color.g + color.b;
                output << "(" << color.r << "," << color.g << "," << color.b << ")" << brightness;
            }
            else { output << "NO_COLOR"; }

            output << ", id=" << id << endl;
            ostringstream retstream;
            retstream << id;
            return retstream.str();
        }
        else
        {
            output << "-- NO_BEACON --" << endl;
            return "";
        }
    }
    catch (NotImplemented const& e)
    {
        output << endl << "NotImplemented while printing info for beacon " << id << ": " << e.what() << endl;
        cerr << endl << "NotImplemented while printing info for beacon " << id << ": " << e.what() << endl;
        return "";
    }
}

MainProgram::CmdResult MainProgram::cmd_find_beacons(ostream& /*output*/, MatchIter begin, MatchIter end)
{
    string name = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto result = ds_.find_beacons(name);
    if (result.empty()) { result = {NO_BEACON}; }

    return {ResultType::IDLIST, result};
}

void MainProgram::test_find_beacons(Stopwatch& watch)
{
    // Choose random number to remove
    if (random_beacons_added_ > 0) // Don't find if there's nothing to find
    {
        auto name = n_to_name(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.find_beacons(name);
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_route_any(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string fromxstr = *begin++;
    string fromystr = *begin++;
    string toxstr = *begin++;
    string toystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int fromx = convert_string_to<int>(fromxstr);
    int fromy = convert_string_to<int>(fromystr);
    int tox = convert_string_to<int>(toxstr);
    int toy = convert_string_to<int>(toystr);

    auto result = ds_.route_any({fromx, fromy}, {tox, toy});

    if (result.empty())
    {
        output << "No path found!" << endl;
    }

    return {ResultType::PATH, result};
}

void MainProgram::test_route_any(Stopwatch& watch)
{
    if (random_beacons_added_ > 0)
    {
        // Choose two random beacons
        auto id1 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        auto id2 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.route_any(ds_.get_coordinates(id1), ds_.get_coordinates(id2));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_route_fastest(ostream& output, MatchIter begin, MatchIter end)
{
    string fromxstr = *begin++;
    string fromystr = *begin++;
    string toxstr = *begin++;
    string toystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int fromx = convert_string_to<int>(fromxstr);
    int fromy = convert_string_to<int>(fromystr);
    int tox = convert_string_to<int>(toxstr);
    int toy = convert_string_to<int>(toystr);

    auto result = ds_.route_fastest({fromx, fromy}, {tox, toy});

    if (result.empty())
    {
        output << "No path found!" << endl;
    }

    return {ResultType::PATH, result};
}

void MainProgram::test_route_fastest(Stopwatch& watch)
{
    if (random_beacons_added_ > 0)
    {
        // Choose two random beacons
        auto id1 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        auto id2 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.route_fastest(ds_.get_coordinates(id1), ds_.get_coordinates(id2));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_route_least_xpoints(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string fromxstr = *begin++;
    string fromystr = *begin++;
    string toxstr = *begin++;
    string toystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int fromx = convert_string_to<int>(fromxstr);
    int fromy = convert_string_to<int>(fromystr);
    int tox = convert_string_to<int>(toxstr);
    int toy = convert_string_to<int>(toystr);

    auto result = ds_.route_least_xpoints({fromx, fromy}, {tox, toy});

    if (result.empty())
    {
        output << "No path found!" << endl;
    }

    return {ResultType::PATH, result};
}

void MainProgram::test_route_least_xpoints(Stopwatch& watch)
{
    if (random_beacons_added_ > 0)
    {
        // Choose two random towns
        auto id1 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        auto id2 = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.route_least_xpoints(ds_.get_coordinates(id1), ds_.get_coordinates(id2));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_route_fibre_cycle(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    string fromxstr = *begin++;
    string fromystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int fromx = convert_string_to<int>(fromxstr);
    int fromy = convert_string_to<int>(fromystr);

    auto result = ds_.route_fibre_cycle({fromx, fromy});

    if (result.empty())
    {
        output << "No fibre cycles found." << endl;
    }

    return {ResultType::CYCLE, result};
}

void MainProgram::test_route_fibre_cycle(Stopwatch& watch)
{
    if (random_beacons_added_ > 0)
    {
        // Choose random town
        auto id = n_to_id(random<decltype(random_beacons_added_)>(0, random_beacons_added_));
        watch.start();
        ds_.route_fibre_cycle(ds_.get_coordinates(id));
        watch.stop();
    }
}

MainProgram::CmdResult MainProgram::cmd_clear_fibres(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    ds_.clear_fibres();
    output << "All fibres removed." << endl;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_randseed(std::ostream& output, MatchIter begin, MatchIter end)
{
    string seedstr = *begin++;
    assert(begin == end && "Invalid number of parameters");

    unsigned long int seed = convert_string_to<unsigned long int>(seedstr);

    rand_engine_.seed(seed);
    init_primes();

    output << "Random seed set to " << seed << endl;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_read(std::ostream& output, MatchIter begin, MatchIter end)
{
    string filename = *begin++;
    string silentstr =  *begin++;
    assert( begin == end && "Impossible number of parameters!");

    bool silent = !silentstr.empty();
    ostream* new_output = &output;

    ostringstream dummystr; // Given as output if "silent" is specified, the output is discarded
    if (silent)
    {
        new_output = &dummystr;
    }

    ifstream input(filename);
    if (input)
    {
        output << "** Commands from '" << filename << "'" << endl;
        command_parser(input, *new_output, PromptStyle::NORMAL);
        if (silent) { output << "...(output discarded in silent mode)..." << endl; }
        output << "** End of commands from '" << filename << "'" << endl;
    }
    else
    {
        output << "Cannot open file '" << filename << "'!" << endl;
    }

    return {};
}


MainProgram::CmdResult MainProgram::cmd_testread(std::ostream& output, MatchIter begin, MatchIter end)
{
    string infilename = *begin++;
    string outfilename = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    ifstream input(infilename);
    if (input)
    {
        ifstream expected_output(outfilename);
        if (output)
        {
            stringstream actual_output;
            command_parser(input, actual_output, PromptStyle::NO_NESTING);

            vector<string> actual_lines;
            while (actual_output)
            {
                string line;
                getline(actual_output, line);
                if (!actual_output) { break; }
                actual_lines.push_back(line);
            }

            vector<string> expected_lines;
            while (expected_output)
            {
                string line;
                getline(expected_output, line);
                if (!expected_output) { break; }
                expected_lines.push_back(line);
            }

            string heading_actual = "Actual output";
            unsigned int actual_max_length = heading_actual.length();
            auto actual_max_iter = max_element(actual_lines.cbegin(), actual_lines.cend(),
                                               [](string s1, string s2){ return s1.length() < s2.length(); });
            if (actual_max_iter != actual_lines.cend())
            {
                actual_max_length = actual_max_iter->length();
            }

            string heading_expected = "Expected output";
            unsigned int expected_max_length = heading_expected.length();
            auto expected_max_iter = max_element(expected_lines.cbegin(), expected_lines.cend(),
                                                    [](string s1, string s2){ return s1.length() < s2.length(); });
            if (expected_max_iter != expected_lines.cend())
            {
                expected_max_length = expected_max_iter->length();
            }

            // New diff, algorithm converted from https://github.com/alexdzyoba/diff
            size_t expsize = expected_lines.size(), actsize = actual_lines.size();
            vector<vector<int>> c(expsize + 1, vector<int>(actsize + 1, 0));

            // Create LCS length matrix
            for (size_t i = 1; i <= expsize; ++i)
            {
                for (size_t j = 1; j <= actsize; ++j)
                {
                    if (expected_lines[i - 1] == actual_lines[j - 1])
                    {
                        c[i][j] = 1 + c[i - 1][j - 1];
                    }
                    else
                    {
                        c[i][j] = max(c[i - 1][j], c[i][j - 1]);
                    }
                }
            }

            size_t i = expsize, j = actsize;
            vector<tuple<char, string, string>> result;

            while (i > 0 || j > 0)
            {
                if (i > 0 && j > 0 && expected_lines[i - 1] == actual_lines[j - 1])
                {
                    result.emplace_back(' ', expected_lines[i - 1], actual_lines[j - 1]);
                    --i; --j;
                }
                else if (i > 0 && j > 0 && c[i][j] == c[i - 1][j - 1])
                {
                    // Lines differ but not part of LCS
                    result.emplace_back('*', expected_lines[i - 1], actual_lines[j - 1]);
                    --i; --j;
                }
                else if (j > 0 && (i == 0 || c[i][j - 1] >= c[i - 1][j]))
                {
                    result.emplace_back('+', "", actual_lines[j - 1]);
                    --j;
                }
                else
                {
                    result.emplace_back('-', expected_lines[i - 1], "");
                    --i;
                }
            }

            reverse(result.begin(), result.end());

            bool lines_ok = true;
            output << "Lines beginning with:" << endl;
            output << "  '+' are extra lines in actual output (not found in expected output)" << endl;
            output << "  '-' are missing lines in actual output (found in expected output)" << endl;
            output << "  '*' are different lines in actual output (compared to expected output)" << endl;
            output << "  " << heading_actual << string(actual_max_length - heading_actual.length(), ' ') << " | " << heading_expected << endl;
            output << "--" << string(actual_max_length, '-') << "-|-" << string(expected_max_length, '-') << endl;
            for (auto const& [diff, expected, actual] : result)
            {
                output << diff << ' ' << actual << string(actual_max_length - actual.length(), ' ')
                       << " | " << expected << endl;
                if (diff != ' ') { lines_ok = false; }
            }

            if (lines_ok)
            {
                output << "**No differences in output.**" << endl;
                if (test_status_ == TestStatus::NOT_RUN)
                {
                    test_status_ = TestStatus::NO_DIFFS;
                }
            }
            else
            {
                output << "**Differences found! (Lines beginning with '+', '-', or '*')**" << endl;
                test_status_ = TestStatus::DIFFS_FOUND;
            }
        }
        else
        {
            output << "Cannot open file '" << outfilename << "'!" << endl;
        }
    }
    else
    {
        output << "Cannot open file '" << infilename << "'!" << endl;
    }

    return {};
}

MainProgram::CmdResult MainProgram::cmd_stopwatch(std::ostream& output, MatchIter begin, MatchIter end)
{
    string on = *begin++;
    string off = *begin++;
    string next = *begin++;
    assert(begin == end && "Invalid number of parameters");

    if (!on.empty())
    {
        stopwatch_mode = StopwatchMode::ON;
        output << "Stopwatch: on" << endl;
    }
    else if (!off.empty())
    {
        stopwatch_mode = StopwatchMode::OFF;
        output << "Stopwatch: off" << endl;
    }
    else if (!next.empty())
    {
        stopwatch_mode = StopwatchMode::NEXT;
        output << "Stopwatch: on for the next command" << endl;
    }
    else
    {
        assert(!"Impossible stopwatch mode!");
    }

    return {};
}

std::string MainProgram::print_coord(Coord coord, std::ostream& output, bool nl)
{
    if (coord != NO_COORD)
    {
        output << "(" << coord.x << "," << coord.y << ")";
        ostringstream retstream;
        retstream << "(" << coord.x << "," << coord.y << ")";
        if (nl) { output << endl; }
        return retstream.str();
    }
    else
    {
        output << "(--NO_COORD--)";
        if (nl) { output << endl; }
        return "";
    }
}

string const beaconidx = "([a-zA-Z0-9-]+)";
string const coordx = "\\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*\\)";
string const optcoordx = "\\([[:space:]]*[0-9]+[[:space:]]*,[[:space:]]*[0-9]+[[:space:]]*\\)";
string const rgbx = "\\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*\\)";
string const numx = "([0-9]+)";
string const wsx = "[[:space:]]+";
string const namex = "([ a-zA-Z0-9-]+)";
string const cmdx = "[0-9a-zA-Z_#]+";

vector<MainProgram::CmdInfo> MainProgram::cmds_ =
{
    {"beacon_info", "ID", beaconidx, &MainProgram::cmd_beacon_info, &MainProgram::test_get_functions },
    {"add_beacon", "ID Name (x,y) (r,g,b)", beaconidx+wsx+namex+wsx+coordx+wsx+rgbx, &MainProgram::cmd_add_beacon, nullptr },
    {"random_add", "number_of_beacons_to_add  (minx,miny) (maxx,maxy) (coordinates optional)",
     numx+"(?:"+wsx+coordx+wsx+coordx+")?", &MainProgram::cmd_random_add, &MainProgram::test_random_add },
    {"random_fibres", "max_number_of_fibres_to_add", numx, &MainProgram::cmd_random_fibres, &MainProgram::test_random_fibres },
    {"random_labyrinth", "xsize ysize extra_routes", numx+wsx+numx+wsx+numx, &MainProgram::cmd_random_labyrinth, nullptr },
    {"all_beacons", "", "", &MainProgram::cmd_all_beacons, &MainProgram::test_all_beacons },
    {"all_xpoints", "", "", &MainProgram::cmd_all_xpoints, &MainProgram::test_all_xpoints },
    {"all_fibres", "", "", &MainProgram::cmd_all_fibres, nullptr },
    {"beacon_count", "", "", &MainProgram::cmd_beacon_count, nullptr },
    {"clear_beacons", "", "", &MainProgram::cmd_clear_beacons, nullptr },
    {"sort_alpha", "", "", &MainProgram::NoParBeaconListCmd<&Datastructures::beacons_alphabetically>, &MainProgram::NoParListTestCmd<&Datastructures::beacons_alphabetically> },
    {"sort_brightness", "", "", &MainProgram::NoParBeaconListCmd<&Datastructures::beacons_brightness_increasing>, &MainProgram::NoParListTestCmd<&Datastructures::beacons_brightness_increasing> },
    {"min_brightness", "", "", &MainProgram::NoParBeaconCmd<&Datastructures::min_brightness>, &MainProgram::NoParBeaconTestCmd<&Datastructures::min_brightness> },
    {"max_brightness", "", "", &MainProgram::NoParBeaconCmd<&Datastructures::max_brightness>, &MainProgram::NoParBeaconTestCmd<&Datastructures::max_brightness> },
    {"find_beacons", "name", namex, &MainProgram::cmd_find_beacons, &MainProgram::test_find_beacons },
    {"change_name", "ID newname", beaconidx+wsx+namex, &MainProgram::cmd_change_name, &MainProgram::test_change_name },
    {"add_lightbeam", "SourceID TargetID", beaconidx+wsx+beaconidx, &MainProgram::cmd_add_lightbeam, nullptr },
    {"lightsources", "BeaconID", beaconidx, &MainProgram::cmd_lightsources, &MainProgram::test_lightsources },
    {"add_fibre", "(x1,y1) (x2,y2) cost", coordx+wsx+coordx+wsx+numx, &MainProgram::cmd_add_fibre, nullptr },
    {"remove_fibre", "(x1,y1) (x2,y2)", coordx+wsx+coordx, &MainProgram::cmd_remove_fibre, &MainProgram::test_remove_fibre },
    {"fibres", "(x,y)", coordx, &MainProgram::cmd_fibres, &MainProgram::test_fibres },
    {"clear_fibres", "", "", &MainProgram::cmd_clear_fibres, nullptr },
    {"path_outbeam", "ID", beaconidx, &MainProgram::cmd_path_outbeam, &MainProgram::test_path_outbeam },
    {"path_inbeam_longest", "ID", beaconidx, &MainProgram::cmd_path_inbeam_longest, &MainProgram::test_path_inbeam_longest },
    {"total_color", "ID", beaconidx, &MainProgram::cmd_total_color, &MainProgram::test_total_color },
    {"route_any", "(x1,y1) (x2,y2)", coordx+wsx+coordx, &MainProgram::cmd_route_any, &MainProgram::test_route_any },
    {"route_fastest", "(x1,y1) (x2,y2)", coordx+wsx+coordx, &MainProgram::cmd_route_fastest, &MainProgram::test_route_fastest },
    {"route_least_xpoints", "(x1,y1) (x2,y2)", coordx+wsx+coordx, &MainProgram::cmd_route_least_xpoints, &MainProgram::test_route_least_xpoints },
    {"route_fibre_cycle", "(x1,y1)", coordx, &MainProgram::cmd_route_fibre_cycle, &MainProgram::test_route_fibre_cycle },
    {"quit", "", "", nullptr, nullptr },
    {"help", "", "", &MainProgram::help_command, nullptr },
    {"read", "\"in-filename\" [silent]", "\"([-a-zA-Z0-9 ./:_]+)\"(?:"+wsx+"(silent))?", &MainProgram::cmd_read, nullptr },
    {"testread", "\"in-filename\" \"out-filename\"", "\"([-a-zA-Z0-9 ./:_]+)\""+wsx+"\"([-a-zA-Z0-9 ./:_]+)\"", &MainProgram::cmd_testread, nullptr },
    {"perftest", "all|compulsory|cmd1[;cmd2...][;extra_add] timeout repeat_count n1[;n2...] (parts in [] are optional, alternatives separated by |)",
     "("+cmdx+"(?:;"+cmdx+")*)"+wsx+numx+wsx+numx+wsx+"([0-9]+(?:;[0-9]+)*)", &MainProgram::cmd_perftest, nullptr },
    {"stopwatch", "on|off|next (alternatives separated by |)", "(?:(on)|(off)|(next))", &MainProgram::cmd_stopwatch, nullptr },
    {"random_seed", "new-random-seed-integer", numx, &MainProgram::cmd_randseed, nullptr },
    {"#", "comment text", ".*", &MainProgram::cmd_comment, &MainProgram::test_comment },
    // {"watchtest", "iterations", numx, &MainProgram::cmd_watchtest, nullptr },
};

void MainProgram::init_cmds()
{
    // Sort commands in alphabetical order
    sort(cmds_.begin(), cmds_.end(), [](auto const& l, auto const& r){ return l.cmd < r.cmd; });
}

// MainProgram::CmdResult MainProgram::cmd_watchtest(std::ostream &output, MatchIter begin, MatchIter end)
// {
//     unsigned int iterations = convert_string_to<unsigned int>(*begin++);
//     assert(begin == end && "Invalid number of parameters");

//     Stopwatch watch(true);
//     for (auto i=0u; i<iterations; ++i)
//     {
//         watch.start();
//         watch.stop();
//     }
//     output << "Stopwatch overhead after " << iterations << " empty iterations: " << endl
//            << "    " << watch.elapsed() << " s, " << watch.elapsed()/iterations << " s/iteration," << endl
//            << "    " << watch.count() << " operations, " << watch.count()/iterations << " operations/iteration" << endl;

//     return {};
// }

MainProgram::CmdResult MainProgram::help_command(std::ostream& output, MatchIter /*begin*/, MatchIter /*end*/)
{
    output << "Commands:" << endl;
    for (auto& i : cmds_)
    {
        output << "  " << i.cmd << " " << i.info << endl;
    }

    return {};
}

MainProgram::CmdResult MainProgram::cmd_perftest(std::ostream& output, MatchIter begin, MatchIter end)
{
#ifdef _GLIBCXX_DEBUG
    output << "WARNING: Debug STL enabled, performance will be worse than expected (maybe also asymptotically)!" << endl;
#endif // _GLIBCXX_DEBUG

    try {
    // Note: everything below is indented too little by one indentation level! (because of try block above)

    vector<string> optional_cmds({"remove_beacon", "path_inbeam_longest", "total_color"});
    vector<string> nondefault_cmds({"#", "all_beacons", "all_xpoints", "remove_beacon", "find_beacons"});

    string commandstr = *begin++;
    unsigned int timeout = convert_string_to<unsigned int>(*begin++);
    unsigned int repeat_count = convert_string_to<unsigned int>(*begin++);
    string sizes = *begin++;
    assert(begin == end && "Invalid number of parameters");

    vector<string> testcmds;
    bool all_cmds = false;
    bool compulsory_cmds = false;
    bool extra_add = false;
    bool additional_get_cmds = false;
    smatch scmd;
    auto cbeg = commandstr.cbegin();
    auto cend = commandstr.cend();
    for ( ; regex_search(cbeg, cend, scmd, commands_regex_); cbeg = scmd.suffix().first)
    {
        testcmds.push_back(scmd[1]);
    }

    if (find(testcmds.begin(), testcmds.end(), "all") != testcmds.end())
    {
        additional_get_cmds = true;
        all_cmds = true;
    }
    if (find(testcmds.begin(), testcmds.end(), "compulsory") != testcmds.end())
    {
        additional_get_cmds = true;
        compulsory_cmds = true;
    }
    auto extra_add_pos = find(testcmds.begin(), testcmds.end(), "extra_add");
    if (extra_add_pos != testcmds.end())
    {
        extra_add = true;
        testcmds.erase(extra_add_pos);
    }

    // Initialize test functions
    vector<void(MainProgram::*)(Stopwatch& watch)> testfuncs;
    if (all_cmds)
    {
        output << "Testing all commands: ";
        for (auto& i : cmds_)
        {
            if (i.testfunc)
            {
                if (find(nondefault_cmds.begin(), nondefault_cmds.end(), i.cmd) == nondefault_cmds.end())
                {
                    output << i.cmd << " ";
                    testfuncs.push_back(i.testfunc);
                }
            }
        }
    }
    else if (compulsory_cmds)
    {
        output << "Testing compulsory commands: ";
        for (auto& i : cmds_)
        {
            if (i.testfunc)
            {
                if (find(nondefault_cmds.begin(), nondefault_cmds.end(), i.cmd) == nondefault_cmds.end() &&
                    find(optional_cmds.begin(), optional_cmds.end(), i.cmd) == optional_cmds.end())
                {
                    output << i.cmd << " ";
                    testfuncs.push_back(i.testfunc);
                }
            }
        }
    }
    else
    {
        output << "Testing commands: ";
        for (auto& i : testcmds)
        {
            auto pos = find_if(cmds_.begin(), cmds_.end(), [&i](auto const& cmd){ return cmd.cmd == i; });
            if (pos != cmds_.end() && pos->testfunc)
            {
                output << i << " ";
                testfuncs.push_back(pos->testfunc);
            }
            else
            {
                output << "(cannot test " << i << ") ";
            }
        }
    }
    output << endl;
    if (extra_add) { output << "Doing 1 random add between commands (not included in measurements)" << endl; }
    output << endl;

    vector<unsigned int> init_ns;
    smatch size;
    auto sbeg = sizes.cbegin();
    auto send = sizes.cend();
    for ( ; regex_search(sbeg, send, size, sizes_regex_); sbeg = size.suffix().first)
    {
        init_ns.push_back(convert_string_to<unsigned int>(size[1]));
    }

    output << "Timeout for each N is " << timeout << " sec. " << endl;
    output << "For each N perform " << repeat_count << " random command(s) from:" << endl;

    if (testfuncs.empty())
    {
        output << "No commands to test!" << endl;
        return {};
    }

#ifdef USE_PERF_EVENT
    output << setw(7) << "N" << ", " << setw(12) << "add (sec)" << ", " << setw(12) << "add (count)" << ", " << setw(12) << "cmds (sec)" << ", "
           << setw(12) << "cmds (count)"  << ", " << setw(12) << "total (sec)" << ", " << setw(12) << "total (count)" << endl;
#else
    output << setw(7) << "N" << ", " << setw(12) << "add (sec)" << ", " << setw(12) << "cmds (sec)" << ", "
           << setw(12) << "total (sec)" << endl;
#endif
    flush_output(output);

    auto stop = false;
    for (unsigned int n : init_ns)
    {
        if (stop) { break; }

        output << setw(7) << n << ", " << flush;

        ds_.clear_beacons();
        ds_.clear_fibres();
        init_primes();

        Stopwatch stopwatch(true); // Stopwatch for actual testing, use also instruction counting, if enabled
        Stopwatch towatch; // Stopwatch for detecting timeouts, no instruction counting
        towatch.start();

        // Add random beacons (+ light sources)
        for (unsigned int i = 0; i < n / 1000; ++i)
        {
            add_random_beacons(1000, &stopwatch);

            if (towatch.elapsed() >= timeout)
            {
                output << "Timeout!" << endl;
                stop = true;
                break;
            }
            if (check_stop())
            {
                output << "Stopped!" << endl;
                stop = true;
                break;
            }
        }
        if (stop) { break; }

        if (n % 1000 != 0)
        {
            add_random_beacons(n % 1000, &stopwatch);
        }

        // Add random fibres
        for (unsigned int i = 0; i < n / 1000; ++i)
        {
            for (unsigned int j=0; j<1000/15; ++j)
            {
                test_random_fibres(stopwatch);
            }

            if (towatch.elapsed() >= timeout)
            {
                output << "Timeout!" << endl;
                stop = true;
                break;
            }
            if (check_stop())
            {
                output << "Stopped!" << endl;
                stop = true;
                break;
            }
        }
        if (stop) { break; }

        if (n % 1000 != 0)
        {
            for (unsigned int j=0; j< (n % 1000)/15; ++j)
            {
                test_random_fibres(stopwatch);
            }
        }

#ifdef USE_PERF_EVENT
        auto addcount = stopwatch.count();
#endif
        auto addsec = stopwatch.elapsed();

#ifdef USE_PERF_EVENT
        output << setw(12) << addsec << ", " << setw(12) << addcount << ", " << flush;
#else
        output << setw(12) << addsec << ", " << flush;
#endif

        for (unsigned int repeat = 0; repeat < repeat_count; ++repeat)
        {
            auto cmdpos = random(testfuncs.begin(), testfuncs.end());

            (this->**cmdpos)(stopwatch);
            if (additional_get_cmds)
            {
                test_get_functions(stopwatch);
            }
            if (extra_add)
            {
                test_extra_add(stopwatch);
            }

            if (repeat % 10 == 0)
            {
                if (towatch.elapsed() >= timeout)
                {
                    output << "Timeout!" << endl;
                    stop = true;
                    break;
                }
                if (check_stop())
                {
                    output << "Stopped!" << endl;
                    stop = true;
                    break;
                }
            }
        }
        if (stop) { break; }

#ifdef USE_PERF_EVENT
        auto totalcount = stopwatch.count();
#endif
        auto totalsec = stopwatch.elapsed();

#ifdef USE_PERF_EVENT
        output << setw(12) << totalsec-addsec << ", " << setw(12) << totalcount-addcount << ", " << setw(12) << totalsec << ", " << setw(12) << totalcount;
#else
        output << setw(12) << totalsec-addsec << ", " << setw(12) << totalsec;
#endif

//        unsigned long int maxmem;
//        string unit;
//        auto [maxmem, unit] = mempeak();
//        maxmem -=  startmem;
//        if (maxmem != 0)
//        {
//            output << ", memory " << maxmem << " " << unit;
//        }
        output << endl;
        flush_output(output);
    }

    output << "Perftest complete." << endl;

    ds_.clear_beacons();
    ds_.clear_fibres();
    init_primes();

    }
    catch (NotImplemented const&)
    {
        // Clean up after NotImplemented
        ds_.clear_beacons();
        ds_.clear_fibres();
        init_primes();
        throw;
    }

#ifdef _GLIBCXX_DEBUG
    output << "WARNING: Debug STL enabled, performance will be worse than expected (maybe also asymptotically)!" << endl;
#endif // _GLIBCXX_DEBUG

    return {};
}

MainProgram::CmdResult MainProgram::cmd_comment(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    return {};
}

// Can be used in the perftest command to skip the actual command (only test adding)
void MainProgram::test_comment(Stopwatch& /*watch*/)
{
}

bool MainProgram::command_parse_line(string inputline, ostream& output)
{
//    static unsigned int nesting_level = 0; // UGLY! Remember nesting level to print correct amount of >:s.
//    if (promptstyle != PromptStyle::NO_NESTING) { ++nesting_level; }

    if (inputline.empty()) { return true; }

    smatch match;
    bool matched = regex_match(inputline, match, cmds_regex_);
    if (matched)
    {
        assert(match.size() == 3);
        string cmd = match[1];
        string params = match[2];

        auto pos = find_if(cmds_.begin(), cmds_.end(), [cmd](CmdInfo const& ci) { return ci.cmd == cmd; });
        assert(pos != cmds_.end());

        smatch match2;
        bool matched2 = regex_match(params, match2, pos->param_regex);
        if (matched2)
        {
            if (pos->func)
            {
                assert(!match2.empty());

                Stopwatch stopwatch;
                bool use_stopwatch = (stopwatch_mode != StopwatchMode::OFF);
                // Reset stopwatch mode if only for the next command
                if (stopwatch_mode == StopwatchMode::NEXT) { stopwatch_mode = StopwatchMode::OFF; }

               TestStatus initial_status = test_status_;
               test_status_ = TestStatus::NOT_RUN;

                if (use_stopwatch)
                {
                    stopwatch.start();
                }

                CmdResult result;
                try
                {
                    result = (this->*(pos->func))(output, ++(match2.begin()), match2.end());
                }
                catch (NotImplemented const& e)
                {
                    output << endl << "NotImplemented from cmd " << pos->cmd << " : " << e.what() << endl;
                    cerr << endl << "NotImplemented from cmd " << pos->cmd << " : " << e.what() << endl;
                }

                if (use_stopwatch)
                {
                    stopwatch.stop();
                }

                switch (result.first)
                {
                    case ResultType::NOTHING:
                    {
                        break;
                    }
                    case ResultType::IDLIST:
                    {
                        auto& res = get<CmdResultIDs>(result.second);
                        if (res.size() == 1 && res.front() == NO_BEACON)
                        {
                            output << "Failed (NO_BEACON returned)!" << endl;
                        }
                        else
                        {
                            if (!res.empty())
                            {
                                if (res.size() == 1) { output << "Beacon:" << endl; }
                                else { output << "Beacons:" << endl; }

                                unsigned int num = 0;
                                for (BeaconID id : res)
                                {
                                    ++num;
                                    if (res.size() > 1) { output << num << ". "; }
                                    else { output << "  "; }
                                    print_beacon(id, output);
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::COORDLIST:
                    {
                        auto& res = get<CmdResultCoords>(result.second);
                        if (res.size() == 1 && res.front() == NO_COORD)
                        {
                            output << "Failed (NO_COORD returned)!" << endl;
                        }
                        else
                        {
                            if (!res.empty())
                            {
                                if (res.size() == 1) { output << "Coord:" << endl; }
                                else { output << "Coords:" << endl; }

                                unsigned int num = 0;
                                for (Coord id : res)
                                {
                                    ++num;
                                    if (res.size() > 1) { output << num << ". "; }
                                    else { output << "  "; }
                                    print_coord(id, output);
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::HIERARCHY:
                    {
                        auto& res = get<CmdResultIDs>(result.second);
                        if (res.size() == 1 && res.front() == NO_BEACON)
                        {
                            output << "Failed (NO_BEACON returned)!" << endl;
                        }
                        else
                        {
                            if (!res.empty())
                            {
                                if (res.size() == 1) { output << "Beacon:" << endl; }
                                else { output << "beacons:" << endl; }

                                unsigned int num = 0;
                                for (BeaconID id : res)
                                {
                                    ++num;
                                    if (res.size() > 1) { output << num << ". "; }
                                    else { output << "  "; }
                                    print_beacon(id, output);
                                }
                            }
                        }
                        break;
                    }
                case ResultType::SOURCELIST:
                {
                    auto res = get<CmdResultIDs>(result.second);
                    assert(res.size() > 0 && "Empty list ind SOURCELIST (not even target id)");
                    auto firstid = res.front();
                    res.erase(res.begin());
                    if (res.size() == 1 && res.front() == NO_BEACON)
                    {
                        output << "Failed (NO_BEACON returned)!" << endl;
                    }
                    else
                    {
                        if (!res.empty())
                        {
                            if (res.size() == 1) { output << "Beacon:" << endl; }
                            else { output << "beacons:" << endl; }

                            unsigned int num = 0;
                            for (BeaconID id : res)
                            {
                                ++num;
                                if (res.size() > 1) { output << num << ". "; }
                                else { output << "  "; }
                                print_beacon(id, output);
                            }
                        }
                    }
                    break;
                }
                    case ResultType::FIBRELIST:
                    {
                        auto& [coord1, res] = get<CmdResultFibres>(result.second);
                        if (!res.empty())
                        {
                            if (res.size() == 1 && get<0>(res.front()) == NO_COORD)
                            {
                                output << "Failed (...NO_COORD... returned)!" << endl;
                            }
                            else
                            {
                                output << "Fibres from ";
                                print_coord(coord1, output, false);
                                output << ":" << endl;
                                unsigned int num = 1;
                                for (auto& r : res)
                                {
                                    auto [coord, cost] = r;
                                    output << num << ". ";
                                    if (coord != NO_COORD)
                                    {
                                        print_coord(coord, output, false);
                                    }
                                    output << " : " << cost;
                                    output << endl;

                                    ++num;
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::PATH:
                    {
                        auto& res = get<CmdResultPath>(result.second);
                        if (!res.empty())
                        {
                            if (res.size() == 1 && get<0>(res.front()) == NO_COORD)
                            {
                                output << "Failed (...NO_COORD... returned)!" << endl;
                            }
                            else
                            {
                                unsigned int num = 1;
                                for (auto& r : res)
                                {
                                    auto [coord, cost] = r;
                                    output << num << ". ";
                                    if (num > 1) { output << "-> "; }
                                    else { output << "   "; }
                                    print_coord(coord, output, false);
                                    output << " : " << cost;
                                    output << endl;

                                    ++num;
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::CYCLE:
                    {
                        auto& res = get<CmdResultCoords>(result.second);
                        if (!res.empty())
                        {
                            if (res.size() == 1 && res.front() == NO_COORD)
                            {
                                output << "Failed (...NO_COORD... returned)!" << endl;
                            }
                            else
                            {
                                unsigned int num = 0;
                                if (res.size() < 2)
                                {
                                    output << "Too small path for cycle!";
                                }
                                else
                                {
                                    auto cycbeg = find(res.begin(), res.end()-1, res.back());
                                    if (cycbeg == res.end())
                                    {
                                        output << "No cycle found in path!";
                                    }
                                    else
                                    {
                                        // Swap cycle so that it starts with smaller id
                                        if (((cycbeg+1) < (res.end()-2)) && (*(res.end()-2) < *(cycbeg+1)))
                                        {
                                            reverse(cycbeg+1, res.end()-1);
                                        }
                                        for (auto i = cycbeg; i != res.end(); ++i)
                                        {
                                            auto coord = *i;
                                            output << num << ". ";
                                            if (num > 0) { output << "-> "; }
                                            else { output << "   "; }
                                            print_coord(coord, output, false);
                                            output << endl;

                                            ++num;
                                        }
                                    }
                                }
                            }
                        }
                        break;
                    }
                    default:
                    {
                        assert(false && "Unsupported result type!");
                    }
                }

                if (result != prev_result)
                {
                    prev_result = std::move(result);
                    view_dirty = true;
                }

                if (use_stopwatch)
                {
                    output << "Command '" << cmd << "': " << stopwatch.elapsed() << " sec" << endl;
                }

                if (test_status_ != TestStatus::NOT_RUN)
                {
                    output << "Testread-tests have been run, " << ((test_status_ == TestStatus::DIFFS_FOUND) ? "differences found!" : "no differences found.") << endl;
                }
                if (test_status_ == TestStatus::NOT_RUN || (test_status_ == TestStatus::NO_DIFFS && initial_status == TestStatus::DIFFS_FOUND))
                {
                    test_status_ = initial_status;
                }
            }
            else
            { // No function to run = quit command
                return false;
            }
        }
        else
        {
            output << "Invalid parameters for command '" << cmd << "'!" << endl;
        }
    }
    else
    {
        output << "Unknown command!" << endl;
    }

    return true; // Signal continuing
}

void MainProgram::command_parser(istream& input, ostream& output, PromptStyle promptstyle)
{
    string line;
    do
    {
//        output << string(nesting_level, '>') << " ";
        output << PROMPT;
        getline(input, line, '\n');

        if (promptstyle != PromptStyle::NO_ECHO)
        {
            output << line << endl;
        }

        if (!input) { break; }

        bool cont = command_parse_line(line, output);
        view_dirty = false; // No need to keep track of individual result changes
        if (!cont) { break; }
    }
    while (input);
    //    if (promptstyle != PromptStyle::NO_NESTING) { --nesting_level; }

    view_dirty = true; // To be safe, assume that results have been changed
}

void MainProgram::setui(MainWindow* ui)
{
    ui_ = ui;
}

#ifdef GRAPHICAL_GUI
void MainProgram::flush_output(std::ostream& output)
{
    if (ui_)
    {
        if (auto soutput = dynamic_cast<ostringstream*>(&output))
        {
            ui_->output_text(*soutput);
        }
    }
}
#else
void MainProgram::flush_output(std::ostream& /*output*/)
{
}
#endif

bool MainProgram::check_stop() const
{
#ifdef GRAPHICAL_GUI
    if (ui_)
    {
        return ui_->check_stop_pressed();
    }
#endif
    return false;
}

std::array<unsigned long int, 20> const MainProgram::primes1{4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,
                                                             5009,   5011,   5021,   5023,   5039,   5051,   5059,   5077,   5081,   5087};
std::array<unsigned long int, 20> const MainProgram::primes2{81031,  81041,  81043,  81047,  81049,  81071,  81077,  81083,  81097,  81101,
                                                             81119,  81131,  81157,  81163,  81173,  81181,  81197,  81199,  81203,  81223};

MainProgram::MainProgram()
{
    rand_engine_.seed(time(nullptr));

    //    startmem = get<0>(mempeak());

    init_primes();
    init_regexs();
    init_cmds();
}

int MainProgram::mainprogram(int argc, char* argv[])
{
    vector<string> args(argv, argv+argc);

    if (args.size() < 1 || args.size() > 3)
    {
        cerr << "Usage: " + ((args.size() > 0) ? args[0] : "<program name>") + " [<command file>|--console|--command command]" << endl;
        return EXIT_FAILURE;
    }

    MainProgram mainprg;

    if (args.size() == 2 && args[1] != "--console")
    {
        string filename = args[1];
        ifstream input(filename);
        if (input)
        {
            mainprg.command_parser(input, cout, MainProgram::PromptStyle::NORMAL);
        }
        else
        {
            cout << "Cannot open file '" << filename << "'!" << endl;
        }
    }
    else if (args.size() == 3 && args[1] == "--command")
    {
        string instr = args[2];
        istringstream input(instr);
        if (input)
        {
            mainprg.command_parser(input, cout, MainProgram::PromptStyle::NORMAL);
        }
        else
        {
            cout << "Cannot use '" << instr << "' as input!" << endl;
        }
    }
    else
    {
        mainprg.command_parser(cin, cout, MainProgram::PromptStyle::NO_ECHO);
    }

    cerr << "Program ended normally." << endl;
    if (mainprg.test_status_ == TestStatus::DIFFS_FOUND)
    {
        return EXIT_FAILURE;
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

void MainProgram::init_primes()
{
    // Initialize id generator
    prime1_ = primes1[random<int>(0, primes1.size())];
    prime2_ = primes2[random<int>(0, primes2.size())];
    random_beacons_added_ = 0;
}

Name MainProgram::n_to_name(unsigned long n)
{
    unsigned long int hash = prime1_*n + prime2_;
    string name;

    while (hash > 0)
    {
        auto hexnum = hash % 26;
        hash /= 26;
        name.push_back('a'+hexnum);
    }

    return name;
}

BeaconID MainProgram::n_to_id(unsigned long n)
{
 ostringstream ostr;
 ostr << "B" << n;
 return ostr.str();
}

Coord MainProgram::n_to_coord(unsigned long n)
{
    unsigned long int hash = prime1_ * n + prime2_;
    hash = hash ^ (hash + 0x9e3779b9 + (hash << 6) + (hash >> 2)); // :-P

    return {static_cast<int>(hash % 1000), static_cast<int>((hash/1000) % 1000)};
}

void MainProgram::create_fibre_labyrinth(std::ostream& output, int xsize, int ysize, int extrafibres)
{
    // Clear existing fibres
    ds_.clear_fibres();

    // Vector to collect all connections
    vector<XpointInfo> xpoints(xsize*ysize);

    // Add random fibre network
    // false = fibre, true = no connection
    vector<bool> fibres(xsize*ysize*ENDDIR);

    vector<long int> tset(xsize*ysize, -1);
    auto fibrenum = fibres.size();
    while (fibrenum > 0)
    {
        // Find a random fibre
        auto fibreiter =  random(fibres.begin(), fibres.end());
        fibreiter = find(fibreiter, fibres.end(), false); // Find next fibre
        if (fibreiter == fibres.end())
        {
            fibreiter = find(fibres.begin(), fibres.end(), false); // Find from beginning
        }
        assert(fibreiter != fibres.end() && "Run out of fibres!");
        *fibreiter = true; // Remove fibre

        // Calculate fibre coordinates and direction
        auto fibreidx = fibreiter - fibres.begin();
        auto [x1, rest] = div(fibreidx, static_cast<decltype(fibreidx)>(ysize*ENDDIR));
        auto [y1, dir] = div(rest, static_cast<decltype(fibreidx)>(ENDDIR));

        // Find the next coordinate
        auto [x2, y2] = move_to_dir({static_cast<int>(x1), static_cast<int>(y1)}, static_cast<Dir>(dir));
        if (x2 < 0 || y2 < 0 || x2 >= xsize || y2 >= ysize) { continue; } // If out of labyrinth area

        auto i1 = x1*ysize + y1;
        auto i2 = x2*ysize + y2;

        auto s1 = i1;
        while (tset[s1] >= 0) { s1 = tset[s1]; }
        auto s2 = i2;
        while (tset[s2] >= 0) { s2 = tset[s2]; }
        if (s1 == s2) { continue; } // Already connected
        tset[s2] += tset[s1];
        tset[s1] = s2;

        xpoints[x1*ysize + y1].dirs[dir] = true;
        xpoints[x2*ysize + y2].dirs[invert_dir[dir]] = true;
        if (tset[s2] == -static_cast<long int>(xsize*ysize)) { break; }
    }

    // Add extra connections
    for (int i=0; i<extrafibres; ++i)
    {
        auto x1 = random(0, xsize);
        auto y1 = random(0, ysize);
        auto dir = random(0, static_cast<int>(ENDDIR));
        auto [x2, y2] = move_to_dir({x1, y1}, static_cast<Dir>(dir));
        if (x2 < 0 || x2 >= xsize || y2 < 0 || y2 >= ysize) { continue; }
        xpoints[x1*ysize + y1].dirs[dir] = true; // New fibre
        xpoints[x2*ysize + y2].dirs[invert_dir[dir]] = true; // And also in other direction
    }

    // Add fibres
    Coord start{random(0, xsize), random(0, ysize)};
    vector<Coord> beaconxys;
    beaconxys.reserve(ds_.beacon_count());
    for (auto id : ds_.all_beacons()) // Make sure all beacons have a corresponding xpoint
    {
        auto [x, y] = ds_.get_coordinates(id);
        auto x2 = x/2;
        auto y2 = y/2;
        if (y % 2 == 0 && x % 2 == y2 % 1)
        {
            beaconxys.push_back({x2, y2});
        }
    }
    sort(beaconxys.begin(), beaconxys.end());
    add_labyrinth_fibres(output, {xsize, ysize}, start, start, ENDDIR, 0, xpoints, beaconxys);
}


Coord MainProgram::move_to_dir(Coord coord, MainProgram::Dir dir)
{
    auto [x2, y2] = coord;
    switch (dir)
    {
        case WEST: { --x2; break; }
        case EAST: { ++x2; break; }
        case NORTHWEST: { ++y2; x2 -= (y2 % 2); break; }
        case NORTHEAST: { x2 += (y2 % 2); ++y2; break; }
        case SOUTHWEST: { --y2; x2 -= (y2 % 2); break; }
        case SOUTHEAST: { x2 += (y2 % 2); --y2; break; }
        default: { assert(!"END as move direction!"); }
    }
    return {x2, y2};
}

void MainProgram::add_labyrinth_fibres(std::ostream& output, Coord size, Coord pos, Coord from, Dir fromdir, Cost fromcost, vector<XpointInfo>& xpoints, vector<Coord> const& beaconxys)
{
    Cost const UNITCOST = 2;
    auto& xpointinfo = xpoints[pos.x*size.y + pos.y];
    if (fromdir!=ENDDIR && xpointinfo.dirs[fromdir] && xpointinfo.dirs.count() == 2 &&
            !binary_search(beaconxys.begin(), beaconxys.end(), pos)) // Only connection is to the same direction, and we are not at a beacon
    {
        xpointinfo.visited = true; // Visited now
        add_labyrinth_fibres(output, size, move_to_dir(pos, fromdir), from, fromdir, fromcost+UNITCOST, xpoints, beaconxys); // Continue in the same direction
    }
    else
    {
        Coord xypos{2*pos.x + (pos.y % 2), 2*pos.y};
        Coord xyfrom{2*from.x + (from.y % 2), 2*from.y};
        if (xyfrom < xypos)
        {
            ds_.add_fibre(xyfrom, xypos, fromcost+UNITCOST); // Add fibre to here
            output << "add_fibre (" << xyfrom.x << "," << xyfrom.y << ") (" << xypos.x << "," << xypos.y << ") " << fromcost+UNITCOST << endl;
        }
        if (xpointinfo.visited) { return; } // Already visited, no need to scan further
        xpointinfo.visited = true; // Visited now
        for (int dir=FIRSTDIR; dir!=ENDDIR; ++dir)
        {
            if (xpointinfo.dirs[dir])
            {
                auto pos2 = move_to_dir(pos, static_cast<Dir>(dir));
                add_labyrinth_fibres(output, size, pos2, pos, static_cast<Dir>(dir), 0, xpoints, beaconxys);
            }
        }
    }
}

void MainProgram::add_random_fibres(std::ostream& output, unsigned int random_fibres)
{
    vector<pair<Coord, Coord>> addedfibres;
    auto beacons = ds_.all_beacons();
    sort(beacons.begin(), beacons.end()); // Sort beacon IDs to get deterministic results

    // Add given number of totally random fibres
    for ( ; random_fibres != 0; --random_fibres)
    {
        auto i1 = random(beacons.begin(), beacons.end());
        auto i2 = random(beacons.begin(), beacons.end());
        auto cost = random(0, 100);
        if (i1 != i2)
        {
            Coord p1 = ds_.get_coordinates(*i1);
            Coord p2 = ds_.get_coordinates(*i2);
            bool intersects = false;
            for (auto const& fibre : addedfibres)
            {
                if (doIntersect(p1, p2, fibre.first, fibre.second)) { intersects = true; break; }
            }
            if (intersects) { continue; }
            auto xyfrom = ds_.get_coordinates(*i1);
            auto xyto = ds_.get_coordinates(*i2);
            ds_.add_fibre(xyfrom, xyto, cost);
            output << "add_fibre (" << xyfrom.x << "," << xyfrom.y << ") (" << xyto.x << "," << xyto.y << ") " << cost << endl;
            addedfibres.push_back({p1, p2});
        }
    }
}

// The functions below are taken and modified from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// point q lies on line segment 'pr'
bool MainProgram::onSegment(Coord p, Coord q, Coord r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
       return true;

    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int MainProgram::orientation(Coord p, Coord q, Coord r)
{
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // colinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool MainProgram::doIntersect(Coord p1, Coord q1, Coord p2, Coord q2)
{
    if ((p1 == p2 && q1 == q2) || (p1 == q2 && q1 == p2)) { return true; } // Same line
    if (p1 == p2 || p1 == q2 || q1 == p2 || q1 == q2) { return false; } // One same point, cannot intersect
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

void MainProgram::init_regexs()
{
    // Create regex <whitespace>(cmd1|cmd2|...)<whitespace>(.*)
    string cmds_regex_str = "[[:space:]]*(";
    bool first = true;
    for (auto& cmd : cmds_)
    {
        cmds_regex_str += (first ? "" : "|") + cmd.cmd;
        first = false;

        cmd.param_regex = regex(cmd.param_regex_str+"[[:space:]]*", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    }
    cmds_regex_str += ")(?:[[:space:]]*$|"+wsx+"(.*))";
    cmds_regex_ = regex(cmds_regex_str, std::regex_constants::ECMAScript | std::regex_constants::optimize);
    coords_regex_ = regex(coordx+"[[:space:]]?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    commands_regex_ = regex("("+cmdx+");?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    sizes_regex_ = regex(numx+";?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
}
