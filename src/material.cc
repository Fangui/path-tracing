#include "material.hh"
#include <iostream>

void Material::dump()
{
  std::cout << "Ns : " << ns << std::endl;
  std::cout << "Ka : " << ka.x_ << " " <<  ka.y_ <<  " " << ka.z_ << std::endl;
  std::cout << "Kd : " << kd.x_ <<  " " << kd.y_ <<  " " << kd.z_ << std::endl;
  std::cout << "Ks : " << ks.x_ <<  " " << ks.y_ <<  " " << ks.z_ << std::endl;
  std::cout << "Ke : " << ke.x_ <<  " " << ke.y_ <<  " " << ke.z_ << std::endl;
  std::cout << "Ni : " << ni << std::endl;
  std::cout << "d : " << d << std::endl;
  std::cout << "illum : " << illum << std::endl;
  std::cout << "Tf : " << tf.x_ << " " <<  tf.y_ <<  " " << tf.z_ << std::endl;
}
