/*
 *       Filename:  doubleVectorConverter.cpp
 *        Created:  09/22/2015 10:44:02 PM
 *         Author:  Zhihao Lou
 */


#include <iostream>
#include <vector>
#include <stdexcept>
#include <unistd.h>
#include <libgen.h>

#include "utils/vectorUtils.h"

int main(int argc, char **argv)
{
    char *inputName = NULL;
    char *outputName = NULL;

    char c;
    char *binname = basename(argv[0]);

    try {
        while ( (c = getopt(argc, argv, "o:")) != -1) {
            switch(c) {
                case 'o':
                    outputName = optarg;
                    break;
                default:
                    throw std::runtime_error("Unrecognized option");
            }
        }
        if (argc <= optind) {
            throw std::runtime_error("Missing input file");
        }
        inputName = argv[optind];
        
        if (!outputName) {
            throw std::runtime_error("Missing output filename");
        }
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr << "Usage: " << binname << " [-o output] input" << std::endl;
        return -1;
    }
    std::vector<double> data;
    readDoubleVectorFromText(data, inputName);
    writeVector(data, outputName);
    return 0;
}

