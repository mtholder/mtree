#if !defined(__MT_INI_OPTIONS_H__)
#define __MT_INI_OPTIONS_H__
#include <algorithm>
#include <map>
#include <cstring>
#include <string>
#include "mt_instance.h"
#include "mt_model_description.h"
namespace mt {
struct CaseIMapComp {
    // slow but portable case insensitive < operator. Not locale sensitive!
    bool operator()(const std::string & lhs, const std::string & rhs) const {
        const auto m = std::min(lhs.length(), rhs.length());
        for (auto i = 0U; i < m; ++i) {
            const int lc = std::tolower(lhs[i]);
            const int rc = std::tolower(rhs[i]);
            if (lc == rc) {
                continue;
            }
            return (lc < rc);
        }
        return lhs.length() < rhs.length();
    }
};
inline std::string compose_illegal_ini_value_error(const std::string & section,
                                                   const std::string & setting,
                                                   const std::string & value) {
    std::string v = "Illegal value \"";
    v.append(value);
    v.append("\" for setting \"");
    v.append(setting);
    v.append("\" in section [");
    v.append(section);
    v.append("].");
    return v;
}
class IllegalINIValueError: public std::runtime_error {
    public:
        IllegalINIValueError(const std::string & section,
                        const std::string & setting,
                        const std::string & value)
            :runtime_error(compose_illegal_ini_value_error(section, setting, value)) {
        }
};

class INIBasedSettings {
    public:
        ProcessActionsEnum action = SCORE_ACTION;
        ModelDescription::AscBiasMode modelAsc = ModelDescription::NO_ASC_BIAS;
};
class INIValueChecker {
    public:
        INIValueChecker() {
            actionAction2Enum["LScore"] = SCORE_ACTION;
            actionAction2Enum["OptimizeBranchLengths"] = OPTIMIZE_BR_LEN;
            actionAction2Enum["OptimizeModel"] = OPTIMIZE_PARS;
            actionAction2Enum["FullOptimization"] = FULL_OPTIMIZE;
            //
            modelAscertainment2Enum["None"] = ModelDescription::NO_ASC_BIAS;
            modelAscertainment2Enum["VarOnly"] = ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS;
            modelAscertainment2Enum["VarOnlyMissing"] = ModelDescription::VAR_ONLY_MISSING_ASC_BIAS;
            modelAscertainment2Enum["ParsInfOnly"] = ModelDescription::PARS_ONLY_NO_MISSING_ASC_BIAS;
            modelAscertainment2Enum["ParsInfOnlyMissing"] = ModelDescription::PARS_ONLY_MISSING_ASC_BIAS;
        }

        bool isLegalActionAction(const std::string & value) const {
            return actionAction2Enum.find(value) != actionAction2Enum.end();
        }
        ProcessActionsEnum parseActionAction(const std::string &value) const {
            auto it = actionAction2Enum.find(value);
            if (it == actionAction2Enum.end()) {
                throw IllegalINIValueError("action", "action", value);
            }
            return it->second;
        }
        ModelDescription::AscBiasMode parseModelAscertainment(const std::string & value) const {
            auto it = modelAscertainment2Enum.find(value);
            if (it == modelAscertainment2Enum.end()) {
                throw IllegalINIValueError("model", "ascertainment", value);
            }
            return it->second;
        }
    private:
        std::map<std::string, ProcessActionsEnum, CaseIMapComp> actionAction2Enum;
        std::map<std::string, ModelDescription::AscBiasMode, CaseIMapComp> modelAscertainment2Enum;
};
class INIReader;
bool iniSettingsAreLegal(INIReader &, const INIValueChecker &, std::ostream &err);

} // namespace
#endif
