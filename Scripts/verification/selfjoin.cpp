#include <atomic>
#include "Dahlia/Tools/Threading.h"
#include <chrono>
#include <cstdio>
int main() {
    DH::Thread t;
    t.Start([&t]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));  // like a real calculation
        t.MarkAsFinished();                                          // exactly what CalculationLayer.cpp:406 does
        std::puts("worker returned normally");
    });
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::puts(t.HasFinished() ? "main: finished" : "main: not finished");
}
