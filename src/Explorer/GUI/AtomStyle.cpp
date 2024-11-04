#include "AtomStyle.h"

namespace SCESS {
    std::unordered_map<std::string_view, AtomStyle> AtomStyle::sAtomStyle =
    {
        {"H", {{1.0f, 1.0f, 1.0f}, 1.0f}},
        {"C", {{0.2f, 0.2f, 0.2f}, 1.5f}},
        {"Li", {{0.8f, 0.9f, 1.0f}, 3.15f}},
        {"N", {{0.0f, 0.5f, 0.5f}, 1.8f}},
        {"O", {{1.0f, 0.0f, 0.0f}, 1.1f}},
        {"B", {{0.100f, 0.714f, 0.757f}, 2.4f}}
    };
}
