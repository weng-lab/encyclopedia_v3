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

        zentlib::BigWig bw(inFnp_);
        for(const auto& chrAndSize : bw.ChromsAndSizes()){
            auto chr = chrAndSize.first;
            bwData[chr] = bw.Data(chr);
        }

        std::vector<double> values;
        for(const auto& chrAndData : bwData){
            for(const auto& d : chrAndData.second){
                for (uint32_t i = d.start; i < d.end; ++i) {
                    if(0){
                        // check unmappable....
                    }
                    values.push_back(d.val);
                }
            }
        }

        // z-score normalize
        a::vec av(values.data(), values.size(), false, true);
        const double mean = a::mean(a::vectorise(av));
        const double stddev = a::stddev(a::vectorise(av));

        for(auto& chrAndData : bwData){
            for(auto& d : chrAndData.second){
                d.val = (d.val - mean) / stddev;
                d.val = clamp(d.val, -1.96, 1.96); // ignore extremes
                d.val += -1.96; // shift up
            }
        }

        bfs::path outFnp = inFnp_.replace_extension(".norm.bed");
        std::ofstream out(outFnp.string());
        if(!out.is_open()){
            throw std::runtime_error("could not open file " +
                                     outFnp.string());
        }
        for(auto& chrAndData : bwData){
            for(auto& d : chrAndData.second){
                out << chrAndData.first << "\t"
                    << d.start << "\t"
                    << d.end << "\t"
                    << d.val << "\n";
            }
        }
        out.close();
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
