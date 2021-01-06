#ifndef __Material_Model_Initializer_h__
#define __Material_Model_Initializer_h__

namespace MatModel
{
	class MaterialModelInitializer
	{
	public:
		static int init();

	private:
		MaterialModelInitializer() {}
		~MaterialModelInitializer();

		static bool is_init;
		static MaterialModelInitializer instance;
	};
}

#endif