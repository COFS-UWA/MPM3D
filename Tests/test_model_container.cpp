#include "Tests_pcp.h"

#include "ModelContainer.h"

using namespace::MatModel;

void test_model_container()
{
	ModelContainer mc;

	mc.add_LinearElasticity(5);
	mc.add_ModifiedCamClay(3);
	mc.add_LinearElasticity(3);
	mc.add_ModifiedCamClay(1);
	mc.add_ModifiedCamClay(1);

	std::cout << "le: " << mc.get_num_LinearElasticity() << "\n";
	std::cout << "mcc: " << mc.get_num_ModifiedCamClay() << "\n";

	double dd = 0.0;
	for (LinearElasticity *le_iter = mc.first_LinearElasticity();
		 mc.is_not_end(le_iter); le_iter = mc.next_LinearElasticity(le_iter))
	{
		le_iter->ext_data_d = dd;
		dd += 1.0;
	}

	for (LinearElasticity *le_iter = mc.first_LinearElasticity();
		mc.is_not_end(le_iter); le_iter = mc.next_LinearElasticity(le_iter))
	{
		std::cout << le_iter->ext_data_d << " ";
	}
	std::cout << "\n";
}
