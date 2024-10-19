
#include "../include/tools/log.h"

namespace ZrRender {
std::string GetFileName(const std::string &) {
    std::string file_path = __FILE__;
    char        cut       = '/';
    file_path.erase(0, file_path.find_last_of(cut) + 1);
    return file_path;
}

}  // namespace ZrRender