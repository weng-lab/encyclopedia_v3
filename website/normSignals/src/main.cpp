#define ARMA_64BIT_WORD
#include <armadillo>

#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <unordered_set>
#include <memory>
#include <unordered_map>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include<json/writer.h>

#include <zi/zargs/zargs.hpp>
ZiARG_string(subPeakFnp, "", "subPeakFnp");
ZiARG_string(bigWigFnp, "", "bigWigFnp");
ZiARG_string(dataFnp, "", "outDir");
ZiARG_bool(debug, false, "debug");
ZiARG_int32(windowSize, 5000, "half of the window size (5000 for histone, 1000 for tf)");
ZiARG_bool(doCorr, false, "do corr");

#include "cpp/files.hpp"
#include "BigWigWrapper.hpp"

#define likely(x) __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)

namespace bib {

namespace a = arma;

class Norm{
    std::string inFnp_;

/*
    void writeSignalOverlapingPeaks(){
        a::mat m(peaks_.size(), distance_, a::fill::zeros);

        // z-score norm
        m = (m - a::mean(a::vectorise(m))) / a::stddev(a::vectorise(m));
        m = a::clamp(m, -1.96, 1.96); // trim

        saveMatrix(m, dataFnp_.string() + ".gz");
    }
*/
public:
    Norm(std::string inFnp)
        : inFnp_(inFnp)
    {}

    void run(){
        if(!bfs::exists(inFnp_)){
            throw std::runtime_error("ERROR: file missing: " + inFnp_);
        }

        std::map<std::string, std::vector<zentlib::Interval>> bwData;

        zentlib::BigWig bw(inFnp_);
        for(const auto& chrAndSize : bw.ChromsAndSizes()){
            auto chr = chrAndSize.first;
            bwData[chr] = bw.Data(chr);

        }


    }
};

} // namespace bib

int main(int argc, char** argv){
    zi::parse_arguments(argc, argv, true);  // modifies argc and argv
    const auto args = std::vector<std::string>(argv + 1, argv + argc);

    bib::Norm n(args.at(0));
    n.run();

    return 0;
}
