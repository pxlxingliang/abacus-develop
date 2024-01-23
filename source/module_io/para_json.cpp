
#include "para_json.h"
#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>


#ifdef __RAPIDJSON
#include "json_output/abacusjson.h"
#endif // __RAPIDJSON

namespace Json
{
void json_output()
{
#ifdef __RAPIDJSON
    gen_general_info();
    AbacusJson::write_to_json("abacus.json");
#endif // __RAPIDJSON
}

void convert_time(std::time_t time_now, std::string& time_str)
{
    std::tm* tm = std::localtime(&time_now);
    std::ostringstream oss;
    oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
    time_str = oss.str();
}

} // namespace Json

