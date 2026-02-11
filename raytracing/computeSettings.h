#pragma once
namespace GOAT
{
	struct computeSettings
	{
		int numThreads = 0; ///< number of threads used for the calculation. Default value is 1, i.e. no parallelization.
	};
}