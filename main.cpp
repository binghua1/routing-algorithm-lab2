#include "ispdData.h"
#include "LayerAssignment.h"
#include "Router2D.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>

namespace {

void preprocessNets(ISPDParser::ispdData &data)
{
    std::vector<ISPDParser::Net *> validNets;
    validNets.reserve(data.nets.size());

    for (ISPDParser::Net *net : data.nets)
    {
        net->pin2D.clear();
        net->pin3D.clear();
        net->twopin.clear();

        for (const auto &pin : net->pins)
        {
            const int x = (std::get<0>(pin) - data.lowerLeftX) / data.tileWidth;
            const int y = (std::get<1>(pin) - data.lowerLeftY) / data.tileHeight;
            const int z = std::get<2>(pin) - 1;

            if (std::none_of(net->pin3D.begin(), net->pin3D.end(), [x, y, z](const ISPDParser::Point &p) {
                    return p.x == x && p.y == y && p.z == z;
                }))
            {
                net->pin3D.emplace_back(x, y, z);
            }

            if (std::none_of(net->pin2D.begin(), net->pin2D.end(), [x, y](const ISPDParser::Point &p) {
                    return p.x == x && p.y == y;
                }))
            {
                net->pin2D.emplace_back(x, y);
            }
        }

        if (net->pin3D.size() > 1000 || net->pin2D.size() <= 1)
        {
            delete net;
            continue;
        }
        validNets.push_back(net);
    }

    data.nets.swap(validNets);
    data.numNet = data.nets.size();
}

} // namespace


int main(int argc, char **argv) {

    assert(argc >= 3 && "Usage: ./router <inputFile> <outputFile>");
    std::ifstream fp(argv[1]);
    assert(fp.is_open() && "Failed to open input file");
    ISPDParser::ispdData *ispdData = ISPDParser::parse(fp);
    fp.close();

    preprocessNets(*ispdData);

    Routing2D::RouterConfig config;
    Routing2D::Router router(*ispdData, config);
    router.build();
    router.route();
    router.exportToIspd();

    LayerAssignment::Graph graph;
    graph.initialLA(*ispdData, 1);
    graph.convertGRtoLA(*ispdData, true);
    graph.COLA(true);
    graph.output3Dresult(argv[2]);

    delete ispdData;
    return 0;
}
