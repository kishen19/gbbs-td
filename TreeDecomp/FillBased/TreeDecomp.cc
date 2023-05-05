// Usage:
// ./TreeDecomp -s com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "TreeDecomp.h"

namespace gbbs {
namespace {

template <class Graph>
double TreeDecomp_runner(Graph& G, commandLine P) {
	int heuristic = P.getOptionIntValue("-h", ALL);
	std::cout << "### Application: Tree Decomposition" << std::endl;
	std::cout << "### Graph: " << P.getArgument(0) << std::endl;
	std::cout << "### Threads: " << num_workers() << std::endl;
	std::cout << "### n: " << G.n << std::endl;
	std::cout << "### m: " << G.m << std::endl;
	std::cout << "### ------------------------------------" << std::endl;
	assert(P.getOption("-s"));

	timer t;
	t.start();
	TreeDecompFill(G, heuristic);
	double tt = t.stop();

	std::cout << "### Running Time: " << tt << std::endl;
	return tt;
}

}  // namespace
}  // namespace gbbs

generate_main(gbbs::TreeDecomp_runner, false);
