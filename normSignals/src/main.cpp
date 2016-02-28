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
#include <json/writer.h>

#include <zi/zargs/zargs.hpp>
ZiARG_string(assembly, "", "assembly");
ZiARG_bool(debug, false, "debug");

#include "cpp/files.hpp"
#include "cpp/utility.hpp"
#include "BigWigWrapper.hpp"

namespace bib {

namespace a = arma;

class Norm{
    bfs::path inFnp_;

public:
    Norm(std::string inFnp)
        : inFnp_(inFnp)
    {}

    void run(){
        if(!bfs::exists(inFnp_)){
            throw std::runtime_error("ERROR: file missing: " +
                                     inFnp_.string());
        }

        std::map<std::string, std::vector<zentlib::Interval>> bwData;

        std::cout << "loading data..." << std::endl;
        zentlib::BigWig bw(inFnp_);
        uint64_t totalBases = 0;
        for(const auto& chrAndSize : bw.ChromsAndSizes()){
            auto chr = chrAndSize.first;
            totalBases += chrAndSize.second + 1;
            bwData[chr] = bw.Data(chr);
        }

        std::cout << "combining all values..." << std::endl;
        std::vector<float> values;
        values.reserve(totalBases);
        for(const auto& chrAndData : bwData){
            auto& intervals = chrAndData.second;
            for(size_t i = 0; i < intervals.size(); ++i){
                const auto& d = intervals[i];
                for (uint32_t i = d.start; i < d.end; ++i) {
                    if(0){
                        // check unmappable....
                    }
                    values.push_back(d.val);
                }
            }
            std::cout << "\t" << chrAndData.first << " "
                      << values.size() << "\n";
        }

        // z-score normalize
        std::cout << "computing mean/stddev..." << std::endl;
        a::fvec av(values.data(), values.size(), false, true);
        const float mean = a::mean(a::vectorise(av));
        const float stddev = a::stddev(a::vectorise(av));

        std::cout << "transforming..." << std::endl;
        for(auto& chrAndData : bwData){
            auto& intervals = chrAndData.second;
            std::cout << "\t" << chrAndData.first << std::endl;
            #pragma omp parallel for
            for(size_t i = 0; i < intervals.size(); ++i){
                auto& d = intervals[i];
                d.val = (d.val - mean) / stddev;
                d.val = clamp(d.val, -1.96, 1.96); // ignore extremes
                d.val += 1.96; // shift up
            }
        }

        std::cout << "writing transformed bed..." << std::endl;
        bfs::path outFnp = inFnp_.replace_extension(".norm.bed");
        outFnp = str::replace(outFnp.string(), "encode/data/", "encode/norm");
        files::ensureDir(outFnp);
        std::ofstream out(outFnp.string());
        if(!out.is_open()){
            throw std::runtime_error("could not open file " +
                                     outFnp.string());
        }
        for(auto& chrAndData : bwData){
            std::cout << "\t" << chrAndData.first << std::endl;
            for(auto& d : chrAndData.second){
                if(0 == d.val){
                    continue;
                }
                out << chrAndData.first << "\t"
                    << d.start << "\t"
                    << d.end << "\t"
                    << d.val << "\n";
            }
        }
        out.close();
        std::cout << "wrote " << outFnp << std::endl;
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
