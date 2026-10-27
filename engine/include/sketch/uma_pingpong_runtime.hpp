#pragma once


#include <array>
#include <cstddef>

namespace Sketch {

enum class UmaSlotState { writable, submitted, completed };

class UmaPingPongContract {
public:
    static constexpr std::size_t slotCount = 2U;

    auto state(std::size_t slot) const -> UmaSlotState {
        return slot < slotCount ? slots_[slot] : UmaSlotState::submitted;
    }

    auto acquireWritable(std::size_t slot) -> bool {
        if (slot >= slotCount || slots_[slot] != UmaSlotState::completed) { return false; }
        slots_[slot] = UmaSlotState::writable;
        return true;
    }

    auto markSubmitted(std::size_t slot) -> bool {
        if (slot >= slotCount || slots_[slot] != UmaSlotState::writable) { return false; }
        slots_[slot] = UmaSlotState::submitted;
        return true;
    }

    auto markCompleted(std::size_t slot) -> bool {
        if (slot >= slotCount || slots_[slot] != UmaSlotState::submitted) { return false; }
        slots_[slot] = UmaSlotState::completed;
        return true;
    }

    auto canPrepareNext(std::size_t currentSlot, std::size_t nextSlot) const -> bool {
        return currentSlot < slotCount && nextSlot < slotCount && currentSlot != nextSlot &&
               slots_[currentSlot] == UmaSlotState::submitted &&
               slots_[nextSlot] == UmaSlotState::completed;
    }

private:
    std::array<UmaSlotState, slotCount> slots_{UmaSlotState::completed, UmaSlotState::completed};
};

}
