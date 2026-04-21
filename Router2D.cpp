#include "Router2D.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>

namespace Routing2D {

namespace {
constexpr double kInfCost = 1e18;
constexpr double kOwnEdgeCost = 0.25;

inline int manhattan(int x1, int y1, int x2, int y2) {
    return std::abs(x1 - x2) + std::abs(y1 - y2);
}

void appendHorizontal(std::vector<std::pair<int, int>> &path, int targetX) {
    while (path.back().first != targetX) {
        const int step = (targetX > path.back().first) ? 1 : -1;
        path.emplace_back(path.back().first + step, path.back().second);
    }
}

void appendVertical(std::vector<std::pair<int, int>> &path, int targetY) {
    while (path.back().second != targetY) {
        const int step = (targetY > path.back().second) ? 1 : -1;
        path.emplace_back(path.back().first, path.back().second + step);
    }
}

} // namespace

Router::Router(ISPDParser::ispdData &data, RouterConfig config)
    : data_(&data), config_(config), gx_(0), gy_(0), numHEdges_(0), numVEdges_(0), approxSubNets_(0) {}

double Router::nowSeconds() {
    using Clock = std::chrono::steady_clock;
    using Seconds = std::chrono::duration<double>;
    return Seconds(Clock::now().time_since_epoch()).count();
}

inline int Router::hIdx(int x, int y) const {
    return y * (gx_ - 1) + x;
}

inline int Router::vIdx(int x, int y) const {
    return y * gx_ + x;
}

void Router::build() {
    gx_ = data_->numXGrid;
    gy_ = data_->numYGrid;
    numHEdges_ = (gx_ - 1) * gy_;
    numVEdges_ = gx_ * (gy_ - 1);

    hCap_.assign(numHEdges_, 0);
    vCap_.assign(numVEdges_, 0);
    hDem_.assign(numHEdges_, 0);
    vDem_.assign(numVEdges_, 0);
    hHist_.assign(numHEdges_, 1.0);
    vHist_.assign(numVEdges_, 1.0);

    buildCapacity();
    buildNetStates();

    if (approxSubNets_ > 100000) {
        config_.maxIterations = std::min(config_.maxIterations, 18);
        config_.rerouteQuotaBase = std::min(config_.rerouteQuotaBase, 1200);
        config_.runtimeBudgetSec = std::min(config_.runtimeBudgetSec, 220.0);
    } else if (approxSubNets_ > 45000) {
        config_.maxIterations = std::min(config_.maxIterations, 22);
        config_.rerouteQuotaBase = std::min(config_.rerouteQuotaBase, 1800);
        config_.runtimeBudgetSec = std::min(config_.runtimeBudgetSec, 250.0);
    }

    if (config_.verbose) {
        std::printf("[Router2D] grid=%dx%d nets=%zu subnets=%lld\n", gx_, gy_, nets_.size(), approxSubNets_);
    }
}

void Router::buildCapacity() {
    const int numLayers = data_->numLayer;
    std::vector<int> wireSize(numLayers, 1);
    for (int layer = 0; layer < numLayers; ++layer) {
        wireSize[layer] = data_->minimumWidth[layer] + data_->minimumSpacing[layer];
        if (wireSize[layer] <= 0) {
            wireSize[layer] = 1;
        }

        const int hTrack = data_->horizontalCapacity[layer] / wireSize[layer];
        const int vTrack = data_->verticalCapacity[layer] / wireSize[layer];
        for (int i = 0; i < numHEdges_; ++i) {
            hCap_[i] += hTrack;
        }
        for (int i = 0; i < numVEdges_; ++i) {
            vCap_[i] += vTrack;
        }
    }

    for (const ISPDParser::CapacityAdj *adj : data_->capacityAdjs) {
        const int x1 = std::get<0>(adj->grid1);
        const int y1 = std::get<1>(adj->grid1);
        const int z1 = std::get<2>(adj->grid1) - 1;
        const int x2 = std::get<0>(adj->grid2);
        const int y2 = std::get<1>(adj->grid2);

        if (z1 < 0 || z1 >= numLayers) {
            continue;
        }

        const bool horizontal = (x1 != x2);
        const int x = std::min(x1, x2);
        const int y = std::min(y1, y2);
        const int oriCap = horizontal ? data_->horizontalCapacity[z1] : data_->verticalCapacity[z1];
        const int adjCap = adj->reducedCapacityLevel;
        const int oriTrack = oriCap / wireSize[z1];
        const int adjTrack = adjCap / wireSize[z1];
        const int delta = oriTrack - adjTrack;

        if (horizontal) {
            if (x >= 0 && x < gx_ - 1 && y >= 0 && y < gy_) {
                hCap_[hIdx(x, y)] -= delta;
            }
        } else {
            if (x >= 0 && x < gx_ && y >= 0 && y < gy_ - 1) {
                vCap_[vIdx(x, y)] -= delta;
            }
        }
    }

    for (int &cap : hCap_) {
        cap = std::max(0, cap);
    }
    for (int &cap : vCap_) {
        cap = std::max(0, cap);
    }
}

void Router::decomposeNetPrim(NetState &netState, const ISPDParser::Net &net) const {
    const std::vector<ISPDParser::Point> &pins = net.pin2D;
    const int n = static_cast<int>(pins.size());
    if (n <= 1) {
        return;
    }
    if (n == 2) {
        netState.subnets.push_back(SubNet{pins[0].x, pins[0].y, pins[1].x, pins[1].y, {}});
        return;
    }

    std::vector<int> parent(n, -1);
    std::vector<int> bestDist(n, std::numeric_limits<int>::max());
    std::vector<char> inTree(n, 0);

    inTree[0] = 1;
    for (int i = 1; i < n; ++i) {
        bestDist[i] = manhattan(pins[0].x, pins[0].y, pins[i].x, pins[i].y);
        parent[i] = 0;
    }

    for (int added = 1; added < n; ++added) {
        int next = -1;
        int nextDist = std::numeric_limits<int>::max();
        for (int i = 0; i < n; ++i) {
            if (!inTree[i] && bestDist[i] < nextDist) {
                nextDist = bestDist[i];
                next = i;
            }
        }

        if (next < 0 || parent[next] < 0) {
            break;
        }

        inTree[next] = 1;
        const int p = parent[next];
        netState.subnets.push_back(SubNet{pins[p].x, pins[p].y, pins[next].x, pins[next].y, {}});

        for (int i = 0; i < n; ++i) {
            if (inTree[i]) {
                continue;
            }
            const int dist = manhattan(pins[next].x, pins[next].y, pins[i].x, pins[i].y);
            if (dist < bestDist[i] || (dist == bestDist[i] && next < parent[i])) {
                bestDist[i] = dist;
                parent[i] = next;
            }
        }
    }
}

void Router::buildNetStates() {
    nets_.clear();
    nets_.reserve(data_->nets.size());

    approxSubNets_ = 0;
    for (int i = 0; i < static_cast<int>(data_->nets.size()); ++i) {
        const ISPDParser::Net &net = *data_->nets[i];
        if (net.pin2D.size() <= 1) {
            continue;
        }

        NetState netState;
        netState.dataIndex = i;

        int minX = std::numeric_limits<int>::max();
        int maxX = std::numeric_limits<int>::min();
        int minY = std::numeric_limits<int>::max();
        int maxY = std::numeric_limits<int>::min();
        for (const ISPDParser::Point &pin : net.pin2D) {
            minX = std::min(minX, pin.x);
            maxX = std::max(maxX, pin.x);
            minY = std::min(minY, pin.y);
            maxY = std::max(maxY, pin.y);
        }
        netState.bbox = (maxX - minX) + (maxY - minY);

        decomposeNetPrim(netState, net);
        approxSubNets_ += static_cast<long long>(netState.subnets.size());

        nets_.push_back(std::move(netState));
    }
}

bool Router::netOwnsEdge(const NetState &net, bool horizontal, int edgeIndex) const {
    const auto &table = horizontal ? net.ownedHEdges : net.ownedVEdges;
    auto it = table.find(edgeIndex);
    return it != table.end() && it->second > 0;
}

double Router::edgeCost(int x1, int y1, int x2, int y2, int iter, const NetState *net) const {
    const bool horizontal = (y1 == y2);
    const int edgeIndex = horizontal ? hIdx(std::min(x1, x2), y1) : vIdx(x1, std::min(y1, y2));

    if (net != nullptr && netOwnsEdge(*net, horizontal, edgeIndex)) {
        return kOwnEdgeCost;
    }

    const int demand = horizontal ? hDem_[edgeIndex] : vDem_[edgeIndex];
    const int cap = horizontal ? hCap_[edgeIndex] : vCap_[edgeIndex];
    const double hist = horizontal ? hHist_[edgeIndex] : vHist_[edgeIndex];

    if (cap <= 0) {
        return 100000.0 + hist * 100.0;
    }

    const int overflow = std::max(0, demand - cap);
    const double util = static_cast<double>(demand + 1) / static_cast<double>(std::max(1, cap));
    const double logistic = 1.0 + 0.8 / (1.0 + std::exp(-2.0 * static_cast<double>(demand + 1 - cap)));
    const double historyPenalty = 1.0 + 0.18 * hist * (1.0 + 0.01 * iter);
    const double overflowPenalty = (overflow > 0) ? (12.0 * overflow * overflow) : 0.0;
    const double nearCapacityPenalty = (util > 1.0) ? (8.0 * (util - 1.0)) : 0.0;

    return logistic * historyPenalty + overflowPenalty + nearCapacityPenalty;
}

double Router::pathCost(const std::vector<GridPoint> &path, int iter, const NetState *net) const {
    if (path.size() <= 1) {
        return 0.0;
    }

    double total = 0.0;
    int bends = 0;
    for (int i = 0; i + 1 < static_cast<int>(path.size()); ++i) {
        const GridPoint &a = path[i];
        const GridPoint &b = path[i + 1];
        total += edgeCost(a.first, a.second, b.first, b.second, iter, net);

        if (i + 2 < static_cast<int>(path.size())) {
            const GridPoint &c = path[i + 2];
            const bool dir1H = (a.second == b.second);
            const bool dir2H = (b.second == c.second);
            if (dir1H != dir2H) {
                ++bends;
            }
        }
    }

    total += bends * 0.45;
    return total;
}

std::vector<Router::GridPoint> Router::buildPathHV(int sx, int sy, int tx, int ty) const {
    std::vector<GridPoint> path;
    path.emplace_back(sx, sy);
    appendHorizontal(path, tx);
    appendVertical(path, ty);
    return path;
}

std::vector<Router::GridPoint> Router::buildPathVH(int sx, int sy, int tx, int ty) const {
    std::vector<GridPoint> path;
    path.emplace_back(sx, sy);
    appendVertical(path, ty);
    appendHorizontal(path, tx);
    return path;
}

std::vector<Router::GridPoint> Router::buildPathZByX(int sx, int sy, int tx, int ty, int midX) const {
    std::vector<GridPoint> path;
    path.emplace_back(sx, sy);
    appendHorizontal(path, midX);
    appendVertical(path, ty);
    appendHorizontal(path, tx);
    return path;
}

std::vector<Router::GridPoint> Router::buildPathZByY(int sx, int sy, int tx, int ty, int midY) const {
    std::vector<GridPoint> path;
    path.emplace_back(sx, sy);
    appendVertical(path, midY);
    appendHorizontal(path, tx);
    appendVertical(path, ty);
    return path;
}

void Router::routePattern(SubNet &subnet, NetState &net, int iter) {
    const int sx = subnet.x1;
    const int sy = subnet.y1;
    const int tx = subnet.x2;
    const int ty = subnet.y2;

    if (sx == tx && sy == ty) {
        subnet.path = {{sx, sy}};
        return;
    }

    std::vector<std::vector<GridPoint>> candidates;
    candidates.reserve(8);
    candidates.push_back(buildPathHV(sx, sy, tx, ty));
    candidates.push_back(buildPathVH(sx, sy, tx, ty));

    if (sx != tx && sy != ty) {
        const int midX = sx + (tx - sx) / 2;
        const int midY = sy + (ty - sy) / 2;
        candidates.push_back(buildPathZByX(sx, sy, tx, ty, midX));
        candidates.push_back(buildPathZByY(sx, sy, tx, ty, midY));

        const int signX = (tx > sx) ? 1 : -1;
        const int signY = (ty > sy) ? 1 : -1;
        if (std::abs(tx - sx) > 1) {
            candidates.push_back(buildPathZByX(sx, sy, tx, ty, sx + signX));
        }
        if (std::abs(ty - sy) > 1) {
            candidates.push_back(buildPathZByY(sx, sy, tx, ty, sy + signY));
        }
    }

    double bestCost = kInfCost;
    int bestIndex = 0;
    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        const double cost = pathCost(candidates[i], iter, &net);
        if (cost < bestCost) {
            bestCost = cost;
            bestIndex = i;
        }
    }

    subnet.path = std::move(candidates[bestIndex]);
}

bool Router::routeMaze(SubNet &subnet, NetState &net, int iter, int margin, double *cost) {
    const int sx = subnet.x1;
    const int sy = subnet.y1;
    const int tx = subnet.x2;
    const int ty = subnet.y2;

    if (sx == tx && sy == ty) {
        subnet.path = {{sx, sy}};
        if (cost != nullptr) {
            *cost = 0.0;
        }
        return true;
    }

    const int x0 = std::max(0, std::min(sx, tx) - margin);
    const int x1 = std::min(gx_ - 1, std::max(sx, tx) + margin);
    const int y0 = std::max(0, std::min(sy, ty) - margin);
    const int y1 = std::min(gy_ - 1, std::max(sy, ty) + margin);

    const int width = x1 - x0 + 1;
    const int height = y1 - y0 + 1;
    const int size = width * height;

    std::vector<double> dist(size, kInfCost);
    std::vector<int> parent(size, -1);

    auto toId = [=](int x, int y) {
        return (y - y0) * width + (x - x0);
    };
    auto fromId = [=](int id) -> GridPoint {
        return {id % width + x0, id / width + y0};
    };

    struct Node {
        double f;
        double g;
        int id;
        bool operator>(const Node &rhs) const {
            return f > rhs.f;
        }
    };

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    const int startId = toId(sx, sy);
    const int targetId = toId(tx, ty);
    dist[startId] = 0.0;
    pq.push(Node{0.0, 0.0, startId});

    static const int dx[4] = {1, -1, 0, 0};
    static const int dy[4] = {0, 0, 1, -1};

    int expansions = 0;
    const int expansionCap = std::min(config_.maxAStarExpansions, size * 8 + 4000);

    while (!pq.empty()) {
        const Node cur = pq.top();
        pq.pop();

        if (cur.id == targetId) {
            break;
        }
        if (cur.g > dist[cur.id] + 1e-9) {
            continue;
        }

        if (++expansions > expansionCap) {
            break;
        }

        const GridPoint p = fromId(cur.id);
        for (int dir = 0; dir < 4; ++dir) {
            const int nx = p.first + dx[dir];
            const int ny = p.second + dy[dir];
            if (nx < x0 || nx > x1 || ny < y0 || ny > y1) {
                continue;
            }

            const int nid = toId(nx, ny);
            const double stepCost = edgeCost(p.first, p.second, nx, ny, iter, &net);
            const double ng = cur.g + stepCost;
            if (ng >= dist[nid] - 1e-9) {
                continue;
            }

            dist[nid] = ng;
            parent[nid] = cur.id;

            // Small admissible heuristic due to own-edge traversal.
            const double h = kOwnEdgeCost * manhattan(nx, ny, tx, ty);
            pq.push(Node{ng + h, ng, nid});
        }
    }

    if (dist[targetId] >= kInfCost / 10.0) {
        return false;
    }

    std::vector<GridPoint> path;
    for (int cur = targetId; cur != -1; cur = parent[cur]) {
        path.push_back(fromId(cur));
    }
    std::reverse(path.begin(), path.end());

    if (path.empty()) {
        return false;
    }

    subnet.path = std::move(path);
    if (cost != nullptr) {
        *cost = dist[targetId];
    }
    return true;
}

void Router::applyPathDemand(const std::vector<GridPoint> &path, int delta, NetState &net) {
    if (path.size() <= 1) {
        return;
    }

    for (int i = 0; i + 1 < static_cast<int>(path.size()); ++i) {
        const int x1 = path[i].first;
        const int y1 = path[i].second;
        const int x2 = path[i + 1].first;
        const int y2 = path[i + 1].second;

        if (y1 == y2) {
            const int idx = hIdx(std::min(x1, x2), y1);
            auto it = net.ownedHEdges.find(idx);
            if (delta > 0) {
                if (it == net.ownedHEdges.end()) {
                    net.ownedHEdges.emplace(idx, 1);
                    ++hDem_[idx];
                } else {
                    ++it->second;
                }
            } else if (it != net.ownedHEdges.end()) {
                --it->second;
                if (it->second <= 0) {
                    net.ownedHEdges.erase(it);
                    --hDem_[idx];
                }
            }
        } else {
            const int idx = vIdx(x1, std::min(y1, y2));
            auto it = net.ownedVEdges.find(idx);
            if (delta > 0) {
                if (it == net.ownedVEdges.end()) {
                    net.ownedVEdges.emplace(idx, 1);
                    ++vDem_[idx];
                } else {
                    ++it->second;
                }
            } else if (it != net.ownedVEdges.end()) {
                --it->second;
                if (it->second <= 0) {
                    net.ownedVEdges.erase(it);
                    --vDem_[idx];
                }
            }
        }
    }
}

void Router::ripUp(SubNet &subnet, NetState &net) {
    if (!subnet.path.empty()) {
        applyPathDemand(subnet.path, -1, net);
        subnet.path.clear();
    }
}

void Router::commit(SubNet &subnet, NetState &net) {
    if (!subnet.path.empty()) {
        applyPathDemand(subnet.path, 1, net);
    }
}

void Router::runInitialPatternRouting() {
    std::vector<int> order(nets_.size());
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (nets_[a].bbox != nets_[b].bbox) {
            return nets_[a].bbox < nets_[b].bbox;
        }
        return nets_[a].subnets.size() < nets_[b].subnets.size();
    });

    for (int index : order) {
        NetState &net = nets_[index];

        std::sort(net.subnets.begin(), net.subnets.end(), [](const SubNet &a, const SubNet &b) {
            return manhattan(a.x1, a.y1, a.x2, a.y2) < manhattan(b.x1, b.y1, b.x2, b.y2);
        });

        for (SubNet &subnet : net.subnets) {
            routePattern(subnet, net, 1);
            commit(subnet, net);
        }
    }
}

int Router::subnetOverflowScore(const SubNet &subnet) const {
    if (subnet.path.size() <= 1) {
        return 0;
    }

    int score = 0;
    for (int i = 0; i + 1 < static_cast<int>(subnet.path.size()); ++i) {
        const int x1 = subnet.path[i].first;
        const int y1 = subnet.path[i].second;
        const int x2 = subnet.path[i + 1].first;
        const int y2 = subnet.path[i + 1].second;

        if (y1 == y2) {
            const int idx = hIdx(std::min(x1, x2), y1);
            if (hDem_[idx] > hCap_[idx]) {
                score += (hDem_[idx] - hCap_[idx]);
            }
        } else {
            const int idx = vIdx(x1, std::min(y1, y2));
            if (vDem_[idx] > vCap_[idx]) {
                score += (vDem_[idx] - vCap_[idx]);
            }
        }
    }

    return score;
}

std::vector<Router::AffectedSubNet> Router::collectAffectedSubnets() {
    std::vector<AffectedSubNet> affected;

    for (NetState &net : nets_) {
        for (SubNet &subnet : net.subnets) {
            const int score = subnetOverflowScore(subnet);
            if (score <= 0) {
                continue;
            }

            affected.push_back(AffectedSubNet{
                score,
                manhattan(subnet.x1, subnet.y1, subnet.x2, subnet.y2),
                net.bbox,
                &net,
                &subnet,
            });
        }
    }

    std::sort(affected.begin(), affected.end(), [](const AffectedSubNet &a, const AffectedSubNet &b) {
        if (a.overflowScore != b.overflowScore) {
            return a.overflowScore > b.overflowScore;
        }
        if (a.bbox != b.bbox) {
            return a.bbox < b.bbox;
        }
        return a.hpwl < b.hpwl;
    });

    return affected;
}

int Router::rerouteQuota(int affectedCount, int iter, double elapsedSec) const {
    if (affectedCount <= 0) {
        return 0;
    }

    int quota = config_.rerouteQuotaBase + iter * config_.rerouteQuotaGrow;
    quota = std::min(quota, affectedCount);

    const double ratio = elapsedSec / std::max(1.0, config_.runtimeBudgetSec);
    if (ratio > 0.75) {
        quota = std::max(200, quota / 3);
    } else if (ratio > 0.55) {
        quota = std::max(300, quota / 2);
    }

    return std::max(1, quota);
}

int Router::marginForIteration(int iter) const {
    int margin = config_.initialBoxMargin + iter * config_.boxGrowPerIter;
    margin = std::min(margin, config_.maxBoxMargin);
    return std::max(2, margin);
}

void Router::updateHistory() {
    for (int i = 0; i < numHEdges_; ++i) {
        if (hDem_[i] > hCap_[i]) {
            hHist_[i] += 1.0;
        } else {
            hHist_[i] = std::max(1.0, hHist_[i] * 0.995);
        }
    }

    for (int i = 0; i < numVEdges_; ++i) {
        if (vDem_[i] > vCap_[i]) {
            vHist_[i] += 1.0;
        } else {
            vHist_[i] = std::max(1.0, vHist_[i] * 0.995);
        }
    }
}

void Router::runIterativeReroute() {
    const double startSec = nowSeconds();

    int previousTOF = totalOverflow();
    int bestTOF = previousTOF;
    int stall = 0;

    for (int iter = 0; iter < config_.maxIterations; ++iter) {
        const double elapsed = nowSeconds() - startSec;
        if (elapsed >= config_.runtimeBudgetSec) {
            if (config_.verbose) {
                std::printf("[Router2D] stop by runtime budget at iter %d\n", iter);
            }
            break;
        }

        const int curTOF = totalOverflow();
        if (curTOF == 0) {
            if (config_.verbose) {
                std::printf("[Router2D] converged at iter %d (TOF=0)\n", iter);
            }
            break;
        }

        std::vector<AffectedSubNet> affected = collectAffectedSubnets();
        if (affected.empty()) {
            break;
        }

        const int quota = rerouteQuota(static_cast<int>(affected.size()), iter, elapsed);
        const int margin = marginForIteration(iter);

        for (int i = 0; i < quota; ++i) {
            AffectedSubNet &cand = affected[i];
            SubNet &subnet = *cand.subnet;
            NetState &net = *cand.net;

            ripUp(subnet, net);

            routePattern(subnet, net, iter + 2);
            const std::vector<GridPoint> patternPath = subnet.path;
            const double patternCost = pathCost(patternPath, iter + 2, &net);

            double mazeCost = kInfCost;
            const bool mazeOK = routeMaze(subnet, net, iter + 2, margin, &mazeCost);
            if (!mazeOK || mazeCost > patternCost * 1.08) {
                subnet.path = patternPath;
            }

            commit(subnet, net);
        }

        updateHistory();

        const int newTOF = totalOverflow();
        if (newTOF < bestTOF) {
            bestTOF = newTOF;
            stall = 0;
        } else {
            ++stall;
        }

        if (config_.verbose) {
            std::printf("[Router2D] iter=%d tof=%d mof=%d wl=%d affected=%zu quota=%d\n",
                        iter,
                        newTOF,
                        maxOverflow(),
                        totalWireLength(),
                        affected.size(),
                        quota);
        }

        if (stall >= config_.stallIterations) {
            if (config_.verbose) {
                std::printf("[Router2D] stop by stall at iter %d\n", iter);
            }
            break;
        }

        if (newTOF >= previousTOF && iter > 2) {
            // Avoid spending too long on non-improving rounds.
            if (stall >= config_.stallIterations / 2) {
                break;
            }
        }

        previousTOF = newTOF;
    }
}

void Router::route() {
    if (nets_.empty()) {
        return;
    }

    const double begin = nowSeconds();
    runInitialPatternRouting();
    runIterativeReroute();

    if (config_.verbose) {
        std::printf("[Router2D] done tof=%d mof=%d wl=%d time=%.2fs\n",
                    totalOverflow(),
                    maxOverflow(),
                    totalWireLength(),
                    nowSeconds() - begin);
    }
}

int Router::totalOverflow() const {
    int tof = 0;
    for (int i = 0; i < numHEdges_; ++i) {
        tof += std::max(0, hDem_[i] - hCap_[i]);
    }
    for (int i = 0; i < numVEdges_; ++i) {
        tof += std::max(0, vDem_[i] - vCap_[i]);
    }
    return tof;
}

int Router::maxOverflow() const {
    int mof = 0;
    for (int i = 0; i < numHEdges_; ++i) {
        mof = std::max(mof, hDem_[i] - hCap_[i]);
    }
    for (int i = 0; i < numVEdges_; ++i) {
        mof = std::max(mof, vDem_[i] - vCap_[i]);
    }
    return mof;
}

int Router::totalWireLength() const {
    int wl = 0;
    for (const NetState &net : nets_) {
        for (const SubNet &subnet : net.subnets) {
            if (!subnet.path.empty()) {
                wl += static_cast<int>(subnet.path.size()) - 1;
            }
        }
    }
    return wl;
}

void Router::exportToIspd() {
    for (NetState &netState : nets_) {
        ISPDParser::Net *net = data_->nets[netState.dataIndex];
        net->twopin.clear();
        net->twopin.reserve(netState.subnets.size());

        for (SubNet &subnet : netState.subnets) {
            if (subnet.path.size() <= 1) {
                subnet.path = buildPathHV(subnet.x1, subnet.y1, subnet.x2, subnet.y2);
            }

            net->twopin.emplace_back();
            ISPDParser::TwoPin &tp = net->twopin.back();
            tp.from = ISPDParser::Point(subnet.x1, subnet.y1);
            tp.to = ISPDParser::Point(subnet.x2, subnet.y2);
            tp.parNet = net;
            tp.path.clear();

            for (int i = 0; i + 1 < static_cast<int>(subnet.path.size()); ++i) {
                const int x1 = subnet.path[i].first;
                const int y1 = subnet.path[i].second;
                const int x2 = subnet.path[i + 1].first;
                const int y2 = subnet.path[i + 1].second;

                if (y1 == y2) {
                    tp.path.emplace_back(std::min(x1, x2), y1, 0, true);
                } else {
                    tp.path.emplace_back(x1, std::min(y1, y2), 0, false);
                }
            }
        }
    }
}

} // namespace Routing2D
