// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)

#ifndef MAINPROGRAM_HH
#define MAINPROGRAM_HH


#ifdef QT_CORE_LIB
#ifndef NO_GRAPHICAL_GUI
#define GRAPHICAL_GUI
#endif
#endif


#include <string>
#include <random>
#include <regex>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <variant>
#include <bitset>
#include <cassert>

#include "datastructures.hh"

using namespace std::string_literals;
static const auto VERSION = "1.0"s;

class MainWindow; // In case there's UI

class MainProgram
{
public:
    MainProgram();


    class Stopwatch;

    enum class PromptStyle { NORMAL, NO_ECHO, NO_NESTING };
    enum class TestStatus { NOT_RUN, NO_DIFFS, DIFFS_FOUND };

    bool command_parse_line(std::string input, std::ostream& output);
    void command_parser(std::istream& input, std::ostream& output, PromptStyle promptstyle);

    void setui(MainWindow* ui);

    void flush_output(std::ostream& output);
    bool check_stop() const;

    static int mainprogram(int argc, char* argv[]);

private:
    Datastructures ds_;
    MainWindow* ui_ = nullptr;

    static std::string const PROMPT;

    std::minstd_rand rand_engine_;

    static std::array<unsigned long int, 20> const primes1;
    static std::array<unsigned long int, 20> const primes2;
    unsigned long int prime1_ = 0; // Will be initialized to random value from above
    unsigned long int prime2_ = 0; // Will be initialized to random value from above
    unsigned long int random_beacons_added_ = 0; // Counter for random beacons added
    void init_primes();
    Name n_to_name(unsigned long int n);
    std::string n_to_id(unsigned long int n);
    Coord n_to_coord(unsigned long int n);


    enum class StopwatchMode { OFF, ON, NEXT };
    StopwatchMode stopwatch_mode = StopwatchMode::OFF;

    enum class ResultType { NOTHING, IDLIST, COORDLIST, FIBRELIST, SOURCELIST, HIERARCHY, PATH, CYCLE };
    using CmdResultIDs = std::vector<BeaconID>;
    using CmdResultCoords = std::vector<Coord>;
    using CmdResultFibres = std::pair<Coord, std::vector<std::pair<Coord, Cost>>>;
    using CmdResultPath = std::vector<std::pair<Coord, Cost>>;
    using CmdResult = std::pair<ResultType, std::variant<CmdResultIDs, CmdResultCoords, CmdResultFibres, CmdResultPath>>;
    CmdResult prev_result;
    bool view_dirty = true;

    TestStatus test_status_ = TestStatus::NOT_RUN;

    using MatchIter = std::smatch::const_iterator;
    struct CmdInfo
    {
        std::string cmd;
        std::string info;
        std::string param_regex_str;
        CmdResult(MainProgram::*func)(std::ostream& output, MatchIter begin, MatchIter end);
        void(MainProgram::*testfunc)(Stopwatch& watch);
        std::regex param_regex = {};
    };
    static std::vector<CmdInfo> cmds_;
    static void init_cmds();

    // Regex objects and their initialization
    std::regex cmds_regex_;
    std::regex coords_regex_;
    std::regex commands_regex_;
    std::regex sizes_regex_;
    void init_regexs();


    CmdResult help_command(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_beacon(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_beacon_info(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_change_name(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_lightbeam(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_fibre(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_remove_fibre(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_path_outbeam(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_path_inbeam_longest(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_total_color(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_random_add(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_randseed(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_read(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_testread(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_beacon_count(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_all_beacons(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_all_xpoints(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_all_fibres(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_lightsources(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_fibres(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_random_fibres(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_random_labyrinth(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_stopwatch(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_clear_beacons(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_find_beacons(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_perftest(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_comment(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_route_any(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_route_fastest(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_route_least_xpoints(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_route_fibre_cycle(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_clear_fibres(std::ostream& output, MatchIter begin, MatchIter end);
    // CmdResult cmd_watchtest(std::ostream& output, MatchIter begin, MatchIter end);

    void test_get_functions(Stopwatch& watch);
    void test_path_outbeam(Stopwatch& watch);
    void test_path_inbeam_longest(Stopwatch& watch);
    void test_total_color(Stopwatch& watch);
    void test_remove_beacon(Stopwatch& watch);
    void test_remove_fibre(Stopwatch& watch);
    void test_change_name(Stopwatch& watch);
    void test_change_color(Stopwatch& watch);
    void test_find_beacons(Stopwatch& watch);
    void test_random_add(Stopwatch& watch);
    void test_all_beacons(Stopwatch& watch);
    void test_all_xpoints(Stopwatch& watch);
    void test_lightsources(Stopwatch& watch);
    void test_fibres(Stopwatch& watch);
    void test_random_fibres(Stopwatch& watch);
    void test_route_any(Stopwatch& watch);
    void test_route_fastest(Stopwatch& watch);
    void test_route_least_xpoints(Stopwatch& watch);
    void test_route_fibre_cycle(Stopwatch& watch);
    void test_comment(Stopwatch& watch);
    void test_extra_add(Stopwatch& watch);

    void add_random_beacons(unsigned int size, Stopwatch* watchp, Coord min = {1,1}, Coord max = {10000, 10000});
    std::string print_beacon(BeaconID id, std::ostream& output);
    std::string print_coord(Coord coord, std::ostream& output, bool nl = true);

    template <typename Type>
    Type random(Type start, Type end);
    template <typename To>
    static To convert_string_to(std::string from);
    template <typename From>
    static std::string convert_to_string(From from);

    template<BeaconID(Datastructures::*MFUNC)()>
    CmdResult NoParBeaconCmd(std::ostream& output, MatchIter begin, MatchIter end);

    template<std::vector<BeaconID>(Datastructures::*MFUNC)()>
    CmdResult NoParBeaconListCmd(std::ostream& output, MatchIter begin, MatchIter end);

    template<BeaconID(Datastructures::*MFUNC)()>
    void NoParBeaconTestCmd(Stopwatch& watch);

    template<std::vector<BeaconID>(Datastructures::*MFUNC)()>
    void NoParListTestCmd(Stopwatch& watch);

    void create_fibre_labyrinth(std::ostream& output, int xsize, int ysize, int extrafibres);

    enum Dir {FIRSTDIR=0, WEST=0, EAST, NORTHWEST, NORTHEAST, SOUTHWEST, SOUTHEAST, ENDDIR};
    std::array<Dir, ENDDIR> invert_dir{{EAST, WEST, SOUTHEAST, SOUTHWEST, NORTHEAST, NORTHWEST}};
    static Coord move_to_dir(Coord coord, Dir dir);
    struct XpointInfo
    {
        std::bitset<6> dirs;
        bool visited;
    };
    void add_labyrinth_fibres(std::ostream& output, Coord size, Coord pos, Coord from, Dir fromdir, Cost fromcost, std::vector<XpointInfo>& xpoints, const std::vector<Coord>& beaconxys);
    void add_random_fibres(std::ostream& output, unsigned int random_fibres);

    // Helper functions for add_random_fibres
    static bool doIntersect(Coord p1, Coord q1, Coord p2, Coord q2);
    static bool onSegment(Coord p, Coord q, Coord r);
    static int orientation(Coord p, Coord q, Coord r);

    friend class MainWindow;
};

template <typename Type>
Type MainProgram::random(Type start, Type end)
{
    auto range = end-start;
    assert(range != 0 && "random() with zero range!");

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine_);

    return static_cast<Type>(start+num);
}

template <typename To>
To MainProgram::convert_string_to(std::string from)
{
    std::istringstream istr(from);
    To result;
    istr >> std::noskipws >> result;
    if (istr.fail() || !istr.eof())
    {
        throw std::invalid_argument("Cannot convert string to required type");
    }
    return result;
}

#ifdef USE_PERF_EVENT
extern "C"
{
#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <cstring>

static long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                int cpu, int group_fd, unsigned long flags)
{
    int ret;

    ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                   group_fd, flags);
    return ret;
}
}
#endif

class MainProgram::Stopwatch
{
public:
    using Clock = std::chrono::high_resolution_clock;

    Stopwatch(bool use_counter = false) : use_counter_(use_counter)
    {
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            memset(&pe_, 0, sizeof(pe_));
            pe_.type = PERF_TYPE_HARDWARE;
            pe_.size = sizeof(pe_);
            pe_.config = PERF_COUNT_HW_INSTRUCTIONS;
            pe_.disabled = 1;
            pe_.exclude_kernel = 1;
            pe_.exclude_hv = 1;

            fd_ = perf_event_open(&pe_, 0, -1, -1, 0);
            if (fd_ == -1) {
                throw "Couldn't open perf events!";
            }
        }
#endif
        reset();
    }

    ~Stopwatch()
    {
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            close(fd_);
        }
#endif
    }

    void start()
    {
        running_ = true;
        starttime_ = Clock::now();
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_RESET, 0);
            read(fd_, &startcount_, sizeof(startcount_));
            ioctl(fd_, PERF_EVENT_IOC_ENABLE, 0);
        }
#endif
    }

    void stop()
    {
        running_ = false;
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0);
            long long c; read(fd_, &c, sizeof(c));
            counter_ += (c - startcount_);
        }
#endif
        elapsed_ += (Clock::now() - starttime_);
    }

    void reset()
    {
        running_ = false;
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0);
            ioctl(fd_, PERF_EVENT_IOC_RESET, 0);
            counter_ = 0;
        }
#endif
        elapsed_ = elapsed_.zero();
    }

    double elapsed()
    {
        if (!running_)
        {
            return static_cast<std::chrono::duration<double>>(elapsed_).count();
        }
        else
        {
            auto total = elapsed_ + (Clock::now() - starttime_);
            return static_cast<std::chrono::duration<double>>(total).count();
        }
    }

#ifdef USE_PERF_EVENT
    long long count()
    {
        if (use_counter_)
        {
            if (!running_)
            {
                return counter_;
            }
            else
            {
                long long c;
                read(fd_, &c, sizeof(c));
                return counter_ + (c - startcount_);
            }
        }
        else
        {
            assert(!"perf_event not enabled during StopWatch creation!");
        }
    }
#endif

private:
    std::chrono::time_point<Clock> starttime_;
    Clock::duration elapsed_ = Clock::duration::zero();
    bool running_ = false;

    [[maybe_unused]] bool use_counter_;
#ifdef USE_PERF_EVENT
    struct perf_event_attr pe_;
    long long startcount_ = 0;
    long long counter_ = 0;
    int fd_ = 0;
#endif
};


template <typename From>
std::string MainProgram::convert_to_string(From from)
{
    std::ostringstream ostr;
    ostr << from;
    if (ostr.fail())
    {
        throw std::invalid_argument("Cannot convert type to string");
    }
    return ostr.str();
}

template<BeaconID(Datastructures::*MFUNC)()>
MainProgram::CmdResult MainProgram::NoParBeaconCmd(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    auto result = (ds_.*MFUNC)();
    return {ResultType::IDLIST, MainProgram::CmdResultIDs{result}};
}

template<std::vector<BeaconID>(Datastructures::*MFUNC)()>
MainProgram::CmdResult MainProgram::NoParBeaconListCmd(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    auto result = (ds_.*MFUNC)();
    return {ResultType::IDLIST, result};
}

template<BeaconID(Datastructures::*MFUNC)()>
void MainProgram::NoParBeaconTestCmd(Stopwatch& watch)
{
    watch.start();
    (ds_.*MFUNC)();
    watch.stop();
}

template<std::vector<BeaconID>(Datastructures::*MFUNC)()>
void MainProgram::NoParListTestCmd(Stopwatch& watch)
{
    watch.start();
    (ds_.*MFUNC)();
    watch.stop();
}


#endif // MAINPROGRAM_HH
