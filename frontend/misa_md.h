//
// Created by gensh(genshenchu@gmail.com) on 2017/10/19.
//

#ifndef MISA_MD_H
#define MISA_MD_H

#include "kiwi_app.h"
#include "simulation.h"
#include "config_parser.h"

class MISAMD : public kiwi::kiwiApp {

public:

    /**
     * Parse command line argv, and print necessary help information (e.g. run: app --help).
     * @param argc argc from function main().
     * @param argv argv from function main().
     * @return false for interrupting the running of the program after parsing argv.
     */
    bool beforeCreate(int argc, char *argv[]) override;

    /**
     * mpi has been initialed, then parse configure file on master processor and synchronize it to other processor.
     * initial architecture environments(e.g. sunway) if necessary.
     */
    void onCreate() override;

    //create boxes and atoms for later simulation.
    //进行区域分解,创建原子
    bool prepare() override;

    /**
     * run simulation.
     */
    void onStart() override;

    void onFinish() override;

    void beforeDestroy() override;

    void onDestroy() override;

private:
    int argc = 0;
    char **argv;

    std::string configFilePath = "config.yaml"; // configure file path default value.
    simulation *pSimulation;
    double wall_clock_begin = 0.0; // the timestamp when program start running

};

#endif //MISA_MD_H
