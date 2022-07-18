// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <random>
#include "Run/Run.h"

#define PY_SSIZE_T_CLEAN
#include "python2.7/Python.h"
#include <stdio.h>

using namespace std;

int     InitializeAll(Run *);
int     LoadConf(string filename, Run *);
int     createfile();

int main(int argc, char *argv[]) {
    //createfile();
    Run * run = new Run();
    InitializeAll(run);
    run->InitializeCellUniform();
    run->InitializeCellDirectors();
    run->updatePolygonVertices();
    run->dumpConfigurationVtk();

//    for (auto cell : run->cells_) {
//        cell->updateVolume();
//        printf("%f\n", cell->volume_);
//    }
//    run->cells_[200]->updateVolume();
//    printf("%f\n", run->cells_[200]->volume_);
//    run->reconnection_->Lth_ = 0.5;
//    run->reconnection_->I_H(run->edges_[2700], true);
//    run->reconnection_->Lth_ = 2.0;
//    run->reconnection_->H_I(run->polygons_[run->polygons_.size()-1], true);

    run->start();

    return 0;
}

int createfile(){
    //char filename[] = "/scripts/tvm/main.py";

    // FILE* PScriptFile = fopen("main.py", "r");
    // if(PScriptFile){
    //     PyRun_SimpleFile(PScriptFile, "main.py");
    //     fclose(PScriptFile);
    // }
    return 0;
}


int InitializeAll(Run * run) {
    printf("Initialization start ...\n");

    // load initial configuration
    ifstream topofile("sample.topo");
    if (!topofile.is_open()) {
        cout << "Error opening sample topo file" << endl;
        exit(1);
    }

    string buffer;
    string delimiter = " ";
    size_t pos = 0;
    long int tmp_id;
    vector<string> tokens;
    vector<vector<string>> lines;


    while (getline(topofile, buffer))
    {
        pos = buffer.find((char)13);
        if (pos != string::npos) {
            buffer = buffer.substr(0, pos);
        }
        if (buffer.length() == 0) continue;

        tokens.clear();
        while ((pos = buffer.find(delimiter)) != string::npos) {
            string token = buffer.substr(0, pos);
            if (token.length() > 0) {
                tokens.push_back(token);
            }
            buffer.erase(0, pos + delimiter.length());
        }
        if (buffer.length() > 0) {
            tokens.push_back(buffer);
        }
        lines.push_back(tokens);
    }

    bool verticesFlag = false;
    bool edgesFlag = false;
    bool polygonsFlag = false;
    bool cellsFlag = false;
    for (int i = 0; i < lines.size(); i++) {
        tokens = lines[i];
        if (tokens[0] == "vertices") {
            verticesFlag = true;
        } else if (tokens[0] == "edges") {
            verticesFlag = false;
            edgesFlag = true;
        } else if (tokens[0] == "polygons") {
            edgesFlag = false;
            polygonsFlag = true;
        } else if (tokens[0] == "cells") {
            polygonsFlag = false;
            cellsFlag = true;
        } else {
            if (verticesFlag) {
                tmp_id = atol(tokens[0].c_str());
                Vertex * vertex = new Vertex(run, tmp_id);
                for (int j = 1; j < tokens.size(); j++) {
                    vertex->position_[j - 1] = atof(tokens[j].c_str());
                }
                run->vertices_.push_back(vertex);
            }
            if (edgesFlag) {
                tmp_id = atol(tokens[0].c_str());
                Edge * edge = new Edge(run, tmp_id);
                for (int j = 1; j < tokens.size(); j++) {
                    tmp_id = atol(tokens[j].c_str());
                    for (auto vertex : run->vertices_) {
                        if (vertex->id_ == tmp_id) {
                            edge->vertices_.push_back(vertex);
                            break;
                        }
                    }
                }
                run->edges_.push_back(edge);
            }
            if (polygonsFlag) {
                tmp_id = atol(tokens[0].c_str());
                Polygon * polygon = new Polygon(run, tmp_id);
                for (int j = 1; j < tokens.size(); j++) {
                    tmp_id = atol(tokens[j].c_str());
                    for (auto edge : run->edges_) {
                        if (edge->id_ == tmp_id) {
                            polygon->edges_.push_back(edge);
                            break;
                        }
                    }
                }
                run->polygons_.push_back(polygon);
            }
            if (cellsFlag) {
                tmp_id = atol(tokens[0].c_str());
                Cell * cell = new Cell(run, tmp_id);
                for (int j = 1; j < tokens.size(); j++) {
                    tmp_id = atol(tokens[j].c_str());
                    for (auto polygon : run->polygons_) {
                        if (polygon->id_ == tmp_id) {
                            cell->polygons_.push_back(polygon);
                            break;
                        }
                    }
                }
                run->cells_.push_back(cell);
            }
        }
    }

    for (auto vertex : run->vertices_) {
        if (run->count_vertices_ < vertex->id_ + 1) {
            run->count_vertices_ = vertex->id_ + 1;
        }
    }
    for (auto edge : run->edges_) {
        if (run->count_edges_ < edge->id_ + 1) {
            run->count_edges_ = edge->id_ + 1;
        }
    }
    for (auto polygon : run->polygons_) {
        if (run->count_polygons_ < polygon->id_ + 1) {
            run->count_polygons_ = polygon->id_ + 1;
        }
    }
    for (auto cell : run->cells_) {
        if (run->count_cells_ < cell->id_ + 1) {
            run->count_cells_ = cell->id_ + 1;
        }
    }
    cout << "Number of vertices: " << run->vertices_.size() << endl;
    cout << "Maximum vertex ID: " << run->count_vertices_ - 1 << endl;
    cout << "Number of edges: " << run->edges_.size() << endl;
    cout << "Maximum edge ID: " << run->count_edges_ - 1 << endl;
    cout << "Number of polygons: " << run->polygons_.size() << endl;
    cout << "Maximum polygon ID: " << run->count_polygons_ - 1 << endl;
    cout << "Number of cells: " << run->cells_.size() << endl;
    cout << "Maximum cell ID: " << run->count_cells_ - 1 << endl;

    // initialize volume object
    run->volume_ = new Volume(run);
    // initialize interface object
    run->interface_ = new Interface(run);
    // initialize reconnection object
    run->reconnection_ = new Reconnection(run);

    LoadConf("conf", run);

    // update geometry and topology information
    run->updateGeoinfo();
    run->updateVertexCells();
    run->volume_->updatePolygonDirections();

    return 0;
}

int LoadConf(string filename, Run * run) {
    ifstream conf(filename.c_str());

    if (!conf.is_open()) {
        cout << "Error opening conf file" << endl;
        exit(1);
    }

    string buffer;
    string delimiter = " ";
    size_t pos = 0;
    vector<string> tokens;
    vector<vector<string>> lines;

    int time_written = 0;
    int dump_written = 0;
    int log_screen_written = 0;
    int s01_written = 0;
    int s06_written = 0;
    int s07_written = 0;
    int s04_written = 0;
    int s05_written = 0;
    int v01_written = 0;
    int v04_written = 0;
    int v05_written = 0;
    int v06_written = 0;
    int v07_written = 0;    
    int Lth_written = 0;
    int L0s_wirrten = 0;
    int temperature_written = 0;
    int temperature_basement = 0;
    int sigma12_written = 0;
    int sigma14_written = 0;
    int sigma15_written = 0;
    int sigma45_written = 0;
    int sigma46_written = 0;
    int sigma56_written = 0;
    int sigma17_written = 0;
    int sigma47_written = 0;
    int sigma57_written = 0;
    int sigma67_written = 0;

    while (getline(conf, buffer))
    {
        pos = buffer.find((char)13);
        if (pos != string::npos) {
            buffer = buffer.substr(0, pos);
        }
        if ((buffer.length() == 0) || (buffer[0] == '#')) continue;

        tokens.clear();
        while ((pos = buffer.find(delimiter)) != string::npos) {
            string token = buffer.substr(0, pos);
            if (token.length() > 0) {
                tokens.push_back(token);
            }
            buffer.erase(0, pos + delimiter.length());
        }
        if (buffer.length() > 0) {
            tokens.push_back(buffer);
        }
        lines.push_back(tokens);
    }


    for (int i = 0; i < lines.size(); i++) {
        tokens = lines[i];
        if (tokens[0] == "time") {
            if (tokens.size() != 4) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->t_start_ = atof(tokens[1].c_str());
            run->t_end_ = atof(tokens[2].c_str());
            run->dt_ = atof(tokens[3].c_str());
            run->dtr_ = 10*run->dt_;
            time_written = 1;
            cout << "time: " << run->t_start_ << " ~ " << run->t_end_ << " ~ " << run->dt_ << " ~ " << run->dtr_ << endl;
        }
        else if (tokens[0] == "dump") {
            if (tokens.size() != 3) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            if (tokens[1] == "vtk") {
                run->dump_period_ = atof(tokens[2].c_str());
                dump_written = 1;
                cout << "dump: " << tokens[1] << " " << run->dump_period_ << endl;
            }
        }
        else if (tokens[0] == "log") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->log_period_ = atof(tokens[1].c_str());
            log_screen_written = 1;
            cout << "log: " << run->log_period_ << endl;
        }
        else if (tokens[0] == "s01") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->s01_ = atof(tokens[1].c_str());
            s01_written = 1;
            cout << "s01: " << run->s01_ << endl;
        }
        else if (tokens[0] == "s06") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->s06_ = atof(tokens[1].c_str());
            s06_written = 1;
            cout << "s06: " << run->s06_ << endl;
        }
        else if (tokens[0] == "s07") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->s07_ = atof(tokens[1].c_str());
            s07_written = 1;
            cout << "s07: " << run->s07_ << endl;
        }
        else if (tokens[0] == "s04") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->s04_ = atof(tokens[1].c_str());
            s04_written = 1;
            cout << "s04: " << run->s04_ << endl;
        }
        else if (tokens[0] == "s05") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->s05_ = atof(tokens[1].c_str());
            s05_written = 1;
            cout << "s05: " << run->s05_ << endl;
        }
         else if (tokens[0] == "v01") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->v01_ = atof(tokens[1].c_str());
            v01_written = 1;
            cout << "v01: " << run->v01_ << endl;
        }
        else if (tokens[0] == "v06") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->v06_ = atof(tokens[1].c_str());
            v06_written = 1;
            cout << "v06: " << run->v06_ << endl;
        }
        else if (tokens[0] == "v07") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->v07_ = atof(tokens[1].c_str());
            v07_written = 1;
            cout << "v07: " << run->v07_ << endl;
        }
        else if (tokens[0] == "v04") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->v04_ = atof(tokens[1].c_str());
            v04_written = 1;
            cout << "v04: " << run->v04_ << endl;
        }
        else if (tokens[0] == "v05") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->v05_ = atof(tokens[1].c_str());
            v05_written = 1;
            cout << "v05: " << run->v05_ << endl;
        }
        else if (tokens[0] == "l0s") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->l0s_ = atof(tokens[1].c_str());
            L0s_wirrten = 1;
            cout << "l0s: " << run->l0s_ << endl;
        }
        else if (tokens[0] == "Lth") {
            if (tokens.size() != 2 && tokens.size() != 3) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->reconnection_->Lth_ = atof(tokens[1].c_str());
            if (tokens.size() == 3) {
                if (tokens[2] == "verbose") {
                    run->reconnection_->verbose_ = true;
                } else {
                    cerr << "conf file error: ";
                    for (int j = 0; j < tokens.size(); j++) {
                        cerr << tokens[j] << " ";
                    }
                    cerr << endl;
                    exit(1);
                }
            }
            Lth_written = 1;
            cout << "Lth: " << run->reconnection_->Lth_ << " verbose: " << run->reconnection_->verbose_ << endl;
        }
        else if (tokens[0] == "T") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->temperature_ = atof(tokens[1].c_str());
            temperature_written = 1;
            cout << "temperature: " << run->temperature_ << endl;
        }
        else if (tokens[0] == "Tbot") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->temperaturebot_ = atof(tokens[1].c_str());
            temperature_basement = 1;
            cout << "temperaturebot: " << run->temperaturebot_ << endl;
        }
        else if (tokens[0] == "sigma12") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma12_ = atof(tokens[1].c_str());
            sigma12_written = 1;
            cout << "sigma12: " << run->sigma12_ << endl;
        }
        else if (tokens[0] == "sigma14") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma14_ = atof(tokens[1].c_str());
            sigma14_written = 1;
            cout << "sigma14: " << run->sigma14_ << endl;
        }
        else if (tokens[0] == "sigma15") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma15_ = atof(tokens[1].c_str());
            sigma15_written = 1;
            cout << "sigma15: " << run->sigma15_ << endl;
        }
        else if (tokens[0] == "sigma45") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma45_ = atof(tokens[1].c_str());
            sigma45_written = 1;
            cout << "sigma45: " << run->sigma45_ << endl;
        }
        else if (tokens[0] == "sigma46") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma46_ = atof(tokens[1].c_str());
            sigma46_written = 1;
            cout << "sigma46: " << run->sigma46_ << endl;
        }
        else if (tokens[0] == "sigma56") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma56_ = atof(tokens[1].c_str());
            sigma56_written = 1;
            cout << "sigma56: " << run->sigma56_ << endl;
        }
        else if (tokens[0] == "sigma17") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma17_ = atof(tokens[1].c_str());
            sigma17_written = 1;
            cout << "sigma17: " << run->sigma17_ << endl;
        }
        else if (tokens[0] == "sigma47") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma47_ = atof(tokens[1].c_str());
            sigma47_written = 1;
            cout << "sigma47: " << run->sigma47_ << endl;
        }
        else if (tokens[0] == "sigma57") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma57_ = atof(tokens[1].c_str());
            sigma57_written = 1;
            cout << "sigma57: " << run->sigma57_ << endl;
        }
        else if (tokens[0] == "sigma67") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->sigma67_ = atof(tokens[1].c_str());
            sigma67_written = 1;
            cout << "sigma67: " << run->sigma67_ << endl;
        }
        else {
            cerr << "conf file error: ";
            for (int j = 0; j < tokens.size(); j++) {
                cerr << tokens[j] << " ";
            }
            cerr << endl;
            exit(1);
        }
    }



    if (conf.bad() || !conf.eof()) {
        cout << "Error reading file [" << filename << "]" << endl;
        exit(1);
    }

    if (time_written == 0) {
        cout << "conf file error: missing time" << endl;
        exit(1);
    }

    if (dump_written == 0) {
        cout << "conf file error: missing dump" << endl;
        exit(1);
    }

    if (log_screen_written == 0) {
        cout << "conf file error: missing log screen" << endl;
        exit(1);
    }

    if (s01_written == 0) {
        cout << "conf file error: s01" << endl;
        exit(1);
    }

    if (v01_written == 0) {
        cout << "conf file error: v01" << endl;
        exit(1);
    }

    if (Lth_written == 0) {
        cout << "conf file error: Lth" << endl;
        exit(1);
    }

    if (temperature_written == 0) {
        cout << "conf file error: temperature" << endl;
        exit(1);
    }

    conf.close();

    return 0;
}