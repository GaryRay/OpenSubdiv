#ifndef CONTEXT_H
#define CONTEXT_H

#include "../common/stopwatch.h"

class Context {
public:
    void BeginTraverse() {
        _traverse.Start();
    }
    void EndTraverse() {
        _traverse.Stop();
    }

    void BeginIntersect() {
        _intersect.Start();
    };
    void EndIntersect() {
        _intersect.Stop();
    }

    void BeginShade() {
        _shade.Start();
    }
    void EndShade() {
        _shade.Stop();
    }

    void Reset() {
        _traverse.Reset();
        _intersect.Reset();
        _shade.Reset();
    }

    double GetTraverseTime() const {
        return _traverse.GetTotalElapsed();
    }

    double GetIntersectTime() const {
        return _intersect.GetTotalElapsed();
    }

    double GetShadeTime() const {
        return _shade.GetTotalElapsed();
    }

private:
    Stopwatch _traverse;
    Stopwatch _intersect;
    Stopwatch _shade;
};


#endif   // CONTEXT_H
