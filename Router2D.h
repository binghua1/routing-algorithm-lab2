#pragma once

#include "ispdData.h"

#include <unordered_map>
#include <utility>
#include <vector>

namespace Routing2D {

struct RouterConfig {
    int maxIterations = 28;
    int stallIterations = 6;
    int initialBoxMargin = 6;
    int boxGrowPerIter = 2;
    int maxBoxMargin = 64;
    int rerouteQuotaBase = 2500;
    int rerouteQuotaGrow = 120;
    int maxAStarExpansions = 120000;
    double runtimeBudgetSec = 280.0;
    bool verbose = true;
};

class Router {
public:
    explicit Router(ISPDParser::ispdData &data, RouterConfig config = RouterConfig());

    void build();
    void route();
    void exportToIspd();

    int totalOverflow() const;
    int maxOverflow() const;
    int totalWireLength() const;

private:
    using GridPoint = std::pair<int, int>;

    struct SubNet {
        int x1 = 0;
        int y1 = 0;
        int x2 = 0;
        int y2 = 0;
        std::vector<GridPoint> path;
    };

    struct NetState {
        int dataIndex = -1;
        int bbox = 0;
        std::vector<SubNet> subnets;
        std::unordered_map<int, int> ownedHEdges;
        std::unordered_map<int, int> ownedVEdges;
    };

    struct AffectedSubNet {
        int overflowScore = 0;
        int hpwl = 0;
        int bbox = 0;
        NetState *net = nullptr;
        SubNet *subnet = nullptr;
    };

    ISPDParser::ispdData *data_;
    RouterConfig config_;

    int gx_;
    int gy_;
    int numHEdges_;
    int numVEdges_;

    std::vector<int> hCap_;
    std::vector<int> vCap_;
    std::vector<int> hDem_;
    std::vector<int> vDem_;
    std::vector<double> hHist_;
    std::vector<double> vHist_;

    std::vector<NetState> nets_;

    long long approxSubNets_;

    static double nowSeconds();

    inline int hIdx(int x, int y) const;
    inline int vIdx(int x, int y) const;

    void buildCapacity();
    void buildNetStates();
    void decomposeNetPrim(NetState &netState, const ISPDParser::Net &net) const;

    void runInitialPatternRouting();
    void runIterativeReroute();

    std::vector<AffectedSubNet> collectAffectedSubnets();
    int subnetOverflowScore(const SubNet &subnet) const;
    int rerouteQuota(int affectedCount, int iter, double elapsedSec) const;
    int marginForIteration(int iter) const;

    void applyPathDemand(const std::vector<GridPoint> &path, int delta, NetState &net);
    void ripUp(SubNet &subnet, NetState &net);
    void commit(SubNet &subnet, NetState &net);
    void updateHistory();

    bool netOwnsEdge(const NetState &net, bool horizontal, int edgeIndex) const;
    double edgeCost(int x1, int y1, int x2, int y2, int iter, const NetState *net) const;
    double pathCost(const std::vector<GridPoint> &path, int iter, const NetState *net) const;

    std::vector<GridPoint> buildPathHV(int sx, int sy, int tx, int ty) const;
    std::vector<GridPoint> buildPathVH(int sx, int sy, int tx, int ty) const;
    std::vector<GridPoint> buildPathZByX(int sx, int sy, int tx, int ty, int midX) const;
    std::vector<GridPoint> buildPathZByY(int sx, int sy, int tx, int ty, int midY) const;

    void routePattern(SubNet &subnet, NetState &net, int iter);
    bool routeMaze(SubNet &subnet, NetState &net, int iter, int margin, double *cost);
};

} // namespace Routing2D
