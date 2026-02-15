// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)

#include "mainprogram.hh"

// Qt generated main window code

#ifdef GRAPHICAL_GUI
#include <QCoreApplication>
#include <QFileDialog>
#include <QDir>
#include <QFont>
#include <QGraphicsItem>
#include <QColor>
#include <QBrush>
#include <QPen>
#include <QGraphicsItem>
#include <QVariant>

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <iostream>

#include <unordered_set>
using std::unordered_set;

#include <unordered_map>
using std::unordered_map;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <variant>
using std::get;

#include <iostream>
using std::cerr;
using std::endl;

#include <sstream>
using std::ostringstream;

#include <cmath>
using std::atan2;

#include <cassert>

#include "mainwindow.hh"
#include "ui_mainwindow.h"

struct XpointInfo
{
    BeaconID id;
    Coord coord;
};

// Needed to be able to store BeaconID in QVariant (in QGraphicsItem)
//Q_DECLARE_METATYPE(BeaconID)
Q_DECLARE_METATYPE(XpointInfo)

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    mainprg_.setui(this);

    // Execute line
    connect(ui->execute_button, &QPushButton::pressed, this, &MainWindow::execute_line);

    // Line edit
    connect(ui->lineEdit, &QLineEdit::returnPressed, this, &MainWindow::execute_line);

    // File selection
    connect(ui->file_button, &QPushButton::pressed, this, &MainWindow::select_file);

    // Command selection
    for (auto& cmd : mainprg_.cmds_)
    {
        ui->cmd_select->addItem(QString::fromStdString(cmd.cmd));
    }
    connect(ui->cmd_select, &QComboBox::activated, this, &MainWindow::cmd_selected);

    // Number selection
    for (auto i = 0; i <= 20; ++i)
    {
        ui->number_select->addItem(QString("%1").arg(i));
    }
    connect(ui->number_select, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), this, &MainWindow::number_selected);

    // Output box
    QFont monofont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    monofont.setStyleHint(QFont::TypeWriter);
    ui->output->setFont(monofont);

    // Initialize graphics scene & view
    gscene_ = new QGraphicsScene(this);
    ui->graphics_view->setScene(gscene_);
    ui->graphics_view->resetTransform();
    ui->graphics_view->setBackgroundBrush(Qt::black);

    // Selecting graphics items by mouse
    connect(gscene_, &QGraphicsScene::selectionChanged, this, &MainWindow::scene_selection_change);
//    connect(this, &MainProgram::signal_clear_selection, this, &MainProgram::clear_selection);

    // Zoom slider changes graphics view scale
    connect(ui->zoom_plus, &QToolButton::clicked, [this]{ this->ui->graphics_view->scale(1.1, 1.1); });
    connect(ui->zoom_minus, &QToolButton::clicked, [this]{ this->ui->graphics_view->scale(1/1.1, 1/1.1); });
    connect(ui->zoom_1, &QToolButton::clicked, [this]{ this->ui->graphics_view->resetTransform(); });
    connect(ui->zoom_fit, &QToolButton::clicked, this, &MainWindow::fit_view);

    // Changing checkboxes updates view
    connect(ui->lightbeams_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->beamcolor_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->beamroutes_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->fibres_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->beacons_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->beaconnames_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);

    // Unchecking lightbeams checkbox disables beam color and routes checkboxes
    connect(ui->lightbeams_checkbox, &QCheckBox::clicked,
            [this]{ this->ui->beamcolor_checkbox->setEnabled(this->ui->lightbeams_checkbox->isChecked()); });
    connect(ui->lightbeams_checkbox, &QCheckBox::clicked,
            [this]{ this->ui->beamroutes_checkbox->setEnabled(this->ui->lightbeams_checkbox->isChecked()); });

    // Unchecking beacons checkbox disables beacon names checkbox
    connect(ui->beacons_checkbox, &QCheckBox::clicked,
            [this]{ this->ui->beaconnames_checkbox->setEnabled(this->ui->beacons_checkbox->isChecked()); });

    // Changing font scale updates view
    connect(ui->fontscale, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::update_view);

    // Changing point scale updates view
    connect(ui->pointscale, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::update_view);

    // Clear input button
    connect(ui->clear_input_button, &QPushButton::clicked, this, &MainWindow::clear_input_line);

    // Stop button
    connect(ui->stop_button, &QPushButton::clicked, [this](){ this->stop_pressed_ = true; });

    clear_input_line();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update_view()
{
    unordered_set<string> errorset;
    try
    {
        gscene_->clear();
        auto pointscale = ui->pointscale->value();
        auto fontscale = ui->fontscale->value();

        // Collect the result of previous operation
        unordered_map<BeaconID, string> result_beacons;
        map<pair<BeaconID, BeaconID>, string> result_hierarchy; // Should be unordered_map, but pair lacks hash...
        unordered_map<Coord, string, CoordHash> result_coords;
        map<pair<Coord, Coord>, string> result_fibres; // Should be unordered_map, but pair lacks hash...
        switch (mainprg_.prev_result.first)
        {
        case MainProgram::ResultType::IDLIST:
            {
                // Copy the id vectors to the result maps
                auto& prev_result = get<MainProgram::CmdResultIDs>(mainprg_.prev_result.second);
                for (unsigned int i=0; i<prev_result.size(); ++i)
                {
                    result_beacons[prev_result[i]] += MainProgram::convert_to_string(i+1)+". ";
                }
            }
            break;
        case MainProgram::ResultType::SOURCELIST:
            {
                // Copy the id vectors to the result maps
                auto& prev_result = get<MainProgram::CmdResultIDs>(mainprg_.prev_result.second);
                BeaconID target_id = NO_BEACON;
                for (unsigned i=0; i<prev_result.size(); ++i)
                {
                    auto id = prev_result[i];
                    if (i != 0)
                    {
                        result_hierarchy[{id, target_id}] += MainProgram::convert_to_string(i)+". ";
                        result_beacons[id] += MainProgram::convert_to_string(i)+". ";
                    }
                    else { target_id = id; }
                }
            }
            break;
        case MainProgram::ResultType::HIERARCHY:
            {
                // Copy the id vectors to the result maps
                auto& prev_result = get<MainProgram::CmdResultIDs>(mainprg_.prev_result.second);
                BeaconID prev_id = NO_BEACON;
                for (unsigned int i=0; i<prev_result.size(); ++i)
                {
                    auto id = prev_result[i];
                    if (i != 0)
                    {
                        result_hierarchy[{prev_id, id}] += MainProgram::convert_to_string(i)+". ";
                    }
                    result_beacons[id] += MainProgram::convert_to_string(i+1)+". ";
                    prev_id = id;
                }
            }
            break;
        case MainProgram::ResultType::COORDLIST:
            {
                // Copy the id vectors to the result maps
                auto& prev_result = get<MainProgram::CmdResultCoords>(mainprg_.prev_result.second);
                for (unsigned int i=0; i<prev_result.size(); ++i)
                {
                    auto coord = prev_result[i];
                    result_coords[coord] += MainProgram::convert_to_string(i+1)+". ";
                }
            }
            break;
        case MainProgram::ResultType::FIBRELIST:
        {
            auto& [target, fibres] = get<MainProgram::CmdResultFibres>(mainprg_.prev_result.second);
            for (unsigned int i=0; i<fibres.size(); ++i)
            {
                auto [coord, cost] = fibres[i];
                auto str = MainProgram::convert_to_string(i+1);
                result_coords[coord] += str+". ";
                result_fibres[{coord, target}] += str+". ";
            }
        }
        break;
        case MainProgram::ResultType::CYCLE:
        {
            auto& prev_result = get<MainProgram::CmdResultCoords>(mainprg_.prev_result.second);
            Coord prev_coord = NO_COORD;
            for (unsigned int i=0; i<prev_result.size(); ++i)
            {
                auto coord = prev_result[i];
                auto str = MainProgram::convert_to_string(i+1);
                result_fibres[{prev_coord,coord}] += str+". ";
                prev_coord = coord;
            }
        }
        break;
        case MainProgram::ResultType::PATH:
            {
                auto& prev_result = get<MainProgram::CmdResultPath>(mainprg_.prev_result.second);
                Coord prev_coord = NO_COORD;
                for (unsigned int i=0; i<prev_result.size(); ++i)
                {
                    auto [coord, cost] = prev_result[i];
                    auto str = MainProgram::convert_to_string(i+1);
                    result_fibres[{prev_coord,coord}] += str+". ";
                    result_coords[coord] += str+". ";
                    prev_coord = coord;
                }
            }
            break;
        case MainProgram::ResultType::NOTHING:
            break;
        default:
            assert(!"Unhandled result type in update_view()!");
        }

        auto beacons = mainprg_.ds_.all_beacons();
        if (beacons.size() == 1 && beacons.front() == NO_BEACON)
        {
            errorset.insert("all_beacons() returned error {NO_BEACON}");
        }

        if (ui->beacons_checkbox->isChecked())
        {
            auto beacons = mainprg_.ds_.all_beacons();
            if (beacons.size() == 1 && beacons.front() == NO_BEACON)
            {
                errorset.insert("Error from GUI: all_beacons() returned error {NO_BEACON}");
                beacons.clear(); // Clear the beacons so that no more errors are caused by NO_BEACON
            }

            for (auto beaconid : beacons)
            {
                QColor beaconcolor = Qt::gray;
                QColor namecolor = Qt::cyan;
                QColor beaconborder = Qt::gray;
                int beaconzvalue = 1;

                try
                {
                    if (beaconid != NO_BEACON)
                    {
                        auto xy = mainprg_.ds_.get_coordinates(beaconid);
                        auto [x,y] = xy;
                        if (x == NO_VALUE || y == NO_VALUE)
                        {
                            errorset.insert("get_coordinates() returned error NO_COORD/NO_VALUE");
                        }

                        if (x == NO_VALUE || y == NO_VALUE)
                        {
                            x = 0; y = 0;
                            beaconcolor = Qt::magenta;
                            namecolor = Qt::magenta;
                            beaconzvalue = 30;
                        }

                        string prefix;
                        auto res_place = result_beacons.find(beaconid);
                        if (res_place != result_beacons.end())
                        {
                            if (result_beacons.size() > 1) { prefix = res_place->second; }
                            namecolor = Qt::red;
                            beaconborder = Qt::red;
                            beaconzvalue = 2;
                        }

                        if (ui->beacons_checkbox->isChecked())
                        {
                            auto groupitem = gscene_->createItemGroup({});
                            groupitem->setFlag(QGraphicsItem::ItemIsSelectable);
                            groupitem->setData(0, QVariant::fromValue(XpointInfo{beaconid, xy}));

                            QPen placepen(beaconborder);
                            placepen.setWidth(0); // Cosmetic pen
                            auto dotitem = gscene_->addEllipse(-4*pointscale, -4*pointscale, 8*pointscale, 8*pointscale,
                                                               placepen, QBrush(beaconcolor));
                            dotitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                            groupitem->addToGroup(dotitem);
                            //        dotitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                            //        dotitem->setData(0, QVariant::fromValue(town));

                            // Draw names
                            string label = prefix;
                            if (ui->beaconnames_checkbox->isChecked())
                            {
                                try
                                {
                                    auto name = mainprg_.ds_.get_name(beaconid);
                                    if (name == NO_NAME)
                                    {
                                        errorset.insert("get_name() returned error NO_NAME");
                                    }

                                    label += name;
                                }
                                catch (NotImplemented const& e)
                                {
                                    errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
                                    cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
                                }
                            }

                            if (!label.empty())
                            {
                                // Create extra item group to be able to set ItemIgnoresTransformations on the correct level (addSimpleText does not allow
                                // setting initial coordinates in item coordinates
                                auto textgroupitem = gscene_->createItemGroup({});
                                auto textitem = gscene_->addSimpleText(QString::fromStdString(label));
                                auto font = textitem->font();
                                font.setPointSizeF(font.pointSizeF()*fontscale);
                                textitem->setFont(font);
                                textitem->setBrush(QBrush(namecolor));
                                textitem->setPos(-textitem->boundingRect().width()/2, -4*pointscale - textitem->boundingRect().height());
                                textgroupitem->addToGroup(textitem);
                                textgroupitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                                groupitem->addToGroup(textgroupitem);
                            }

                            groupitem->setPos(20*x, -20*y);
                            groupitem->setZValue(beaconzvalue);
                        }

                        // Draw light beams
                        if (ui->lightbeams_checkbox->isChecked())
                        {
                            auto sources = mainprg_.ds_.get_lightsources(beaconid);
                            for (auto& sourceid : sources)
                            {
                                try
                                {
                                    auto [vx,vy] = mainprg_.ds_.get_coordinates(sourceid);
                                    QColor beamcolor = Qt::white;
                                    Qt::PenStyle penstyle = Qt::SolidLine;
                                    int zvalue = -2;

                                    if (vx == NO_VALUE || vy == NO_VALUE)
                                    {
                                        errorset.insert("NO_COORD/NO_VALUE from get_coordinates while updating graphics: ");
                                        cerr << endl << "NO_COORD/NO_VALUE from get_coordinates while updating graphics: " << endl;
                                        continue;
                                    }

                                    if (ui->beamcolor_checkbox->isChecked())
                                    {
                                        auto totcol = mainprg_.ds_.total_color(sourceid);
                                        if (totcol != NO_COLOR)
                                        {
                                            beamcolor = QColor(totcol.r, totcol.g, totcol.b);
                                        }
                                    }
                                    auto linecolor = beamcolor;
                                    auto arrowcolor = beamcolor;

                                    if (result_hierarchy.find({sourceid, beaconid}) != result_hierarchy.end())
                                    {
                                        linecolor = Qt::red;
                                        penstyle = Qt::DotLine;
                                        zvalue = 10;
                                    }

                                    // auto pen = QPen(linecolor);
                                    // pen.setWidth(0); // "Cosmetic" pen
                                    // auto lineitem = gscene_->addLine(20*x, -20*y, 20*rx, -20*ry, pen);
                                    auto pen = QPen(linecolor, 3*pointscale);
                                    auto arrowpen = QPen(arrowcolor);
                                    pen.setStyle(penstyle);

                                    if (ui->beamroutes_checkbox->isChecked())
                                    {
                                        auto route = mainprg_.ds_.route_any({vx, vy}, {x, y});
                                        if (!(route.empty() || route.front().first == NO_COORD))
                                        {
                                            route.pop_back(); // Remove the last segment
                                            for (auto& fibre : route)
                                            {
                                                auto [x2, y2] = fibre.first;
                                                QLineF line(QPointF(20*x2, -20*y2), QPointF(20*vx, -20*vy));
                                                auto lineitem = gscene_->addLine(line, pen);
                                                lineitem->setZValue(zvalue);
                                                vx = x2;
                                                vy = y2;
                                            }
                                        }
                                    }
                                    QLineF line(QPointF(20*x, -20*y), QPointF(20*vx, -20*vy));
                                    auto lineitem = gscene_->addLine(line, pen);
                                    lineitem->setZValue(zvalue);

                                    double const PI  = 3.141592653589793238463;
                                    auto arrowSize = 15*pointscale;

                                    double angle = atan2(-line.dy(), line.dx());

                                    QPointF arrowP1 = line.p1() + QPointF(sin(angle + PI / 3) * arrowSize,
                                                                    cos(angle + PI / 3) * arrowSize);
                                    QPointF arrowP2 = line.p1() + QPointF(sin(angle + PI - PI / 3) * arrowSize,
                                                                    cos(angle + PI - PI / 3) * arrowSize);

                                    QPolygonF arrowHead;
                                    arrowHead << line.p1() << arrowP1 << arrowP2;
                                    auto headitem = gscene_->addPolygon(arrowHead, arrowpen, QBrush(arrowpen.color()));
                                    headitem->setZValue(zvalue);
                                }
                                catch (NotImplemented const& e)
                                {
                                    errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
                                    cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
                                }
                            }
                        }
                    }
                }
                catch (NotImplemented const& e)
                {
                    errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
                    cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
                }
            }
        }

        // Draw fibres
        if (ui->fibres_checkbox->isChecked())
        {
            try
            {
                auto xpoints = mainprg_.ds_.all_xpoints();
                if (xpoints.size() == 1 && xpoints.front() == NO_COORD)
                {
                    errorset.insert("all_xpoints() returned error {NO_COORD}");
                }

                for (auto& xy1 : xpoints)
                {
                    QColor dotcolor = Qt::black;
                    QColor labelcolor = Qt::cyan;
                    QColor dotborder = Qt::gray;
                    int dotzvalue = -1;

                    string label;
                    auto res_place = result_coords.find(xy1);
                    if (res_place != result_coords.end())
                    {
                        if (result_coords.size() > 1) { label = res_place->second; }
                        dotcolor = Qt::red;
                        dotborder = Qt::red;
                        labelcolor = Qt::red;
                        dotzvalue = 2;
                    }

                    auto [x,y] = xy1;
                    if (x == NO_VALUE || y == NO_VALUE)
                    {
                        errorset.insert("NO_COORD/NO_VALUE in xpoints while updating graphics: ");
                        cerr << endl << "NO_COORD/NO_VALUE in xpoints while updating graphics: " << endl;
                        continue;
                    }

                    try
                    {
                        auto groupitem = gscene_->createItemGroup({});
                        groupitem->setFlag(QGraphicsItem::ItemIsSelectable);
                        groupitem->setData(0, QVariant::fromValue(XpointInfo{NO_BEACON, xy1}));

                        QPen xroadpen(dotborder);
                        xroadpen.setWidth(0); // Cosmetic pen
                        auto dotitem2 = gscene_->addEllipse(-5*pointscale, -5*pointscale, 10*pointscale, 10*pointscale,
                                                           xroadpen, QBrush(dotcolor));
                        dotitem2->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                        groupitem->addToGroup(dotitem2);
                        //        dotitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                        //        dotitem->setData(0, QVariant::fromValue(town));

                        // Draw label
                        if (!label.empty())
                        {
                            // Create extra item group to be able to set ItemIgnoresTransformations on the correct level (addSimpleText does not allow
                            // setting initial coordinates in item coordinates
                            auto textgroupitem = gscene_->createItemGroup({});
                            auto textitem = gscene_->addSimpleText(QString::fromStdString(label));
                            auto font = textitem->font();
                            font.setPointSizeF(font.pointSizeF()*fontscale);
                            textitem->setFont(font);
                            textitem->setBrush(QBrush(labelcolor));
                            textitem->setPos(-textitem->boundingRect().width()/2, +4*pointscale /*- textitem->boundingRect().height()*/);
                            textgroupitem->addToGroup(textitem);
                            textgroupitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                            groupitem->addToGroup(textgroupitem);
                        }

                        groupitem->setPos(20*x, -20*y);
                        groupitem->setZValue(dotzvalue);

                        dotitem2->setFlag(QGraphicsItem::ItemIsSelectable);
                        dotitem2->setData(0, QVariant::fromValue(XpointInfo{NO_BEACON, xy1}));

                        auto fibres = mainprg_.ds_.get_fibres_from(xy1);
                        for (auto& fibre : fibres)
                        {
                            auto xy2 = fibre.first;
                            auto [rx,ry] = xy2;

                            if (rx == NO_VALUE || ry == NO_VALUE)
                            {
                                errorset.insert("NO_COORD/NO_VALUE from get_fibres_from while updating graphics: ");
                                cerr << endl << "NO_COORD/NO_VALUE from get_fibres_from while updating graphics: " << endl;
                                continue;
                            }

                            QColor linecolor = Qt::gray;
                            int zvalue = -2;

                            if (result_fibres.find({xy1, xy2}) != result_fibres.end())
                            {
                                linecolor = Qt::red;
                                zvalue = 10;
                            }

                            auto pen = QPen(linecolor);
                            pen.setWidth(0); // "Cosmetic" pen
                            auto lineitem = gscene_->addLine(20*x, -20*y, 20*rx, -20*ry, pen);
                            lineitem->setZValue(zvalue);
                        }
                    }
                    catch (NotImplemented const& e)
                    {
                        errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
                        cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
                    }
                }
            }
            catch (NotImplemented const& e)
            {
                errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
                cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
            }
        }

        gscene_->setSceneRect(gscene_->itemsBoundingRect());
    }
    catch (NotImplemented const& e)
    {
        errorset.insert(string("NotImplemented while updating graphics: ") + e.what());
        cerr << endl << "NotImplemented while updating graphics: " << e.what() << endl;
    }

    if (!errorset.empty())
    {
        ostringstream errorstream;
        for (auto const& errormsg : errorset)
        {
            errorstream << "Error from GUI: " << errormsg << endl;
        }
        output_text(errorstream);
        output_text_end();
    }
}

void MainWindow::output_text(ostringstream& output)
{
    string outstr = output.str();
    if (!outstr.empty())
    {
        if (outstr.back() == '\n') { outstr.pop_back(); } // Remove trailing newline
        ui->output->appendPlainText(QString::fromStdString(outstr));
        ui->output->ensureCursorVisible();
        ui->output->repaint();
    }

    output.str(""); // Clear the stream, because it has already been output
}

void MainWindow::output_text_end()
{
    ui->output->moveCursor(QTextCursor::End);
    ui->output->ensureCursorVisible();
    ui->output->repaint();
}

bool MainWindow::check_stop_pressed() const
{
    QCoreApplication::processEvents();
    return stop_pressed_;
}

void MainWindow::execute_line()
{
    auto line = ui->lineEdit->text();
    clear_input_line();
    ui->output->appendPlainText(QString::fromStdString(MainProgram::PROMPT)+line);

    ui->execute_button->setEnabled(false);
    ui->stop_button->setEnabled(true);
    stop_pressed_ = false;

    ostringstream output;
    bool cont = mainprg_.command_parse_line(line.toStdString(), output);
    ui->lineEdit->clear();
    output_text(output);
    output_text_end();

    ui->stop_button->setEnabled(false);
    ui->execute_button->setEnabled(true);
    stop_pressed_ = false;

    ui->lineEdit->setFocus();

//    if (mainprg_.view_dirty)
//    {
//        update_view();
//        mainprg_.view_dirty = false;
//    }
    update_view();

    if (!cont)
    {
        close();
    }
}

void MainWindow::cmd_selected(int idx)
{
    ui->lineEdit->insert(QString::fromStdString(mainprg_.cmds_[idx].cmd+" "));
    ui->cmd_info_text->setText(QString::fromStdString(mainprg_.cmds_[idx].cmd+" "+mainprg_.cmds_[idx].info));

    ui->lineEdit->setFocus();
}

void MainWindow::number_selected(int idx)
{
    ui->lineEdit->insert(ui->number_select->itemText(idx)+" ");

    ui->lineEdit->setFocus();
}

void MainWindow::select_file()
{
    QFileDialog dialog(this, "Select file", "", "Command list (*.txt)");
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if (dialog.exec())
    {
        auto filenames = dialog.selectedFiles();
        for (auto& i : filenames)
        {
            auto filename = QDir("").relativeFilePath(i);
            ui->lineEdit->insert("\""+filename+"\" ");
        }
    }

    ui->lineEdit->setFocus();
}

void MainWindow::clear_input_line()
{
    ui->cmd_info_text->clear();
    ui->lineEdit->clear();
    ui->lineEdit->setFocus();
}

void MainWindow::fit_view()
{
    ui->graphics_view->fitInView(gscene_->itemsBoundingRect(), Qt::KeepAspectRatio);
}

void MainWindow::scene_selection_change()
{
    auto items = gscene_->selectedItems();
    if (!items.empty())
    {
        for (auto i : items)
        {
            auto [beaconid, coord] = i->data(0).value<XpointInfo>();
//            assert(beaconid!=NO_ID && "NO_ID as beacon in graphics view");
            if (!selection_clear_in_progress)
            {
                ostringstream output;
                output << "*click* ";
                if (ui->pick_beacon->isChecked() && beaconid != NO_BEACON)
                {
                    auto beaconstr = mainprg_.print_beacon(beaconid, output);
                    if (!beaconstr.empty()) { ui->lineEdit->insert(QString::fromStdString(beaconstr+" ")); }
                }
                else
                {
                    auto coordstr = mainprg_.print_coord(coord, output);
                    if (!coordstr.empty()) { ui->lineEdit->insert(QString::fromStdString(coordstr)+" "); }
                }
                output_text(output);
                output_text_end();
            }
            i->setSelected(false);
            selection_clear_in_progress = !selection_clear_in_progress;
        }
    }
}

void MainWindow::clear_selection()
{
    gscene_->clearSelection();
}



// Originally in main.cc
#include <QApplication>

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        // try to find tests-dir
        QDir current = QDir::current();
        auto count = 10;
        do {
            if( current.cd("tests") )
                break;
            current.cd("..");
        } while( --count > 0 );
        auto selected = current.path();
        if( QDir::setCurrent( selected ) )
            qDebug() << "Using tests from:" << selected;



        QApplication a(argc, argv);
        MainWindow w;
        w.show();

        auto status = a.exec();
        cerr << "Program ended normally." << endl;
        return status;
    }
    else
    {
        return MainProgram::mainprogram(argc, argv);
    }
}

#else

int main(int argc, char *argv[])
{
    return MainProgram::mainprogram(argc, argv);
}
#endif
