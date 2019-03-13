#include <memory>

#include "network_coding.hxx"
#include "catch.hpp"

using namespace ff;

TEST_CASE("Tests basic usage of linear combination generator", "[nc_generator_basic]")
{
  //auto rng = std::make_shared<simple_rng>();
  simple_rng rng;

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 5);
  
  nc_generator g(&rng, &datasets, classes);

  std::vector<linear_combination> cs;
  for(size_t i=0; i<10; i++) {
    auto c = g.new_linear_combination();
    cs.push_back(c);
  }

  // decode!
  nc_absorber decoder(&rng, 5);
  for(size_t i=0; i<10; i++) {
    bool fin = decoder.inject(cs.at(i));
    if(i < 4) {
      REQUIRE ( fin == false );
    }
    else {
      REQUIRE ( fin == true );
    }
    // printf("i,fin: %u,%u\n",i,fin);
  }
  auto result = decoder.solve();

  for(size_t i=0; i<datasets.size(); i++) {
    REQUIRE( datasets.at(i) == result.at(i) );
  }
}

TEST_CASE("Tests 100x1000 system for correct decoding", "[nc_100x1000]")
{
  std::mt19937 engine(0);
  std::uniform_int_distribution<short> dist(0,255);

  uint32_t generation_size = 100;
  
  std::shared_ptr<gf_rng> rng = std::make_shared<simple_rng>();
  std::vector<std::vector<uint8_t> > data(generation_size);

  // generate data
  for(size_t i=0; i<generation_size; i++) {
    data.at(i) = std::vector<uint8_t>(1000);
    for(size_t j=0; j<1000; j++) {
      // fill data with random bytes
      data.at(i).at(j) = dist(engine);
    }
  }
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 5);
  	
  nc_generator g(rng.get(), &data, classes);
  nc_absorber a(rng.get(),generation_size);
  
  // generate linear combinations until decode is possible
  size_t i=0;
  bool fin;
  do {
    linear_combination lc = g.new_linear_combination();
    fin = a.inject(lc);
    i++;
  } while(!fin);
  REQUIRE (i == 100);

  auto solution = a.solve();

  for(size_t i=0; i<data.size(); i++) {
    REQUIRE( solution.at(i) == data.at(i) );
  }
}

TEST_CASE("Tests partial decoding", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 5);
  
  nc_generator g(rng.get(), &datasets, classes);
  nc_absorber decoder( rng.get(), 5);

  bool fin;
  REQUIRE ( decoder.get_current_degree() == 0 );
  fin = decoder.inject(g.new_zerotail_linear_combination(1));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 1 );
  fin = decoder.inject(g.new_zerotail_linear_combination(4));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 2 );
  fin = decoder.inject(g.new_zerotail_linear_combination(2));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 3 );
  fin = decoder.inject(g.new_zerotail_linear_combination(5));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 4 );
  fin = decoder.inject(g.new_zerotail_linear_combination(3));
  REQUIRE ( fin == true );
  REQUIRE ( decoder.get_current_degree() == 5 );

  auto result = decoder.solve();

  // check result
  for(size_t i=0; i<datasets.size(); i++) {
    for(size_t j=0; j<datasets.at(0).size(); j++) {
      REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
    }
  }
}

TEST_CASE("Tests partial 3x3 decoding", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 5);
  
  nc_generator g(rng.get(), &datasets, classes);
  nc_absorber decoder(rng.get(), datasets.size());

  bool fin;
  REQUIRE ( decoder.get_current_degree() == 0 );
  fin = decoder.inject(g.new_zerotail_linear_combination(2));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 1 );
  fin = decoder.inject(g.new_zerotail_linear_combination(3));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 2 );
  fin = decoder.inject(g.new_zerotail_linear_combination(1));
  REQUIRE ( fin == true );
  REQUIRE ( decoder.get_current_degree() == 3 );

  auto result = decoder.solve();

  // check result
  for(size_t i=0; i<datasets.size(); i++) {
    for(size_t j=0; j<datasets.at(0).size(); j++) {
      REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
    }
  }
}

TEST_CASE("Tests partial 3x3 decoding (2)", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});

  std::vector<size_t> classes;
  classes.push_back((size_t) 3);
    
  nc_generator g(rng.get(), &datasets, classes);
  nc_absorber decoder(rng.get(), datasets.size());

  bool fin;
  REQUIRE ( decoder.get_current_degree() == 0 );
  fin = decoder.inject(g.new_zerotail_linear_combination(3));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 1 );
  fin = decoder.inject(g.new_zerotail_linear_combination(1));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 2 );
  fin = decoder.inject(g.new_zerotail_linear_combination(2));
  REQUIRE ( fin == true );
  REQUIRE ( decoder.get_current_degree() == 3 );

  auto result = decoder.solve();

  // check result
  for(size_t i=0; i<datasets.size(); i++) {
    for(size_t j=0; j<datasets.at(0).size(); j++) {
      REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
    }
  }
}

TEST_CASE("Tests partial 4x4 decoding", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({87,12,127});

  std::vector<size_t> classes;
  classes.push_back((size_t) 5);
    
  nc_generator g(rng.get(), &datasets, classes);
  nc_absorber decoder(rng.get(), datasets.size());

  bool fin;
  REQUIRE ( decoder.get_current_degree() == 0 );
  fin = decoder.inject(g.new_zerotail_linear_combination(4));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 1 );
  fin = decoder.inject(g.new_zerotail_linear_combination(2));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 2 );
  fin = decoder.inject(g.new_zerotail_linear_combination(3));
  REQUIRE ( fin == false );
  REQUIRE ( decoder.get_current_degree() == 3 );
  fin = decoder.inject(g.new_zerotail_linear_combination(1));
  REQUIRE ( fin == true );
  REQUIRE ( decoder.get_current_degree() == 4 );

  auto result = decoder.solve();

  // check result
  for(size_t i=0; i<datasets.size(); i++) {
    for(size_t j=0; j<datasets.at(0).size(); j++) {
      REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
    }
  }
}

TEST_CASE("Tests various partial decodings and crazy counter elimination", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  for(size_t dim = 5; dim < 20; dim++) {
    printf("d=%zu\n",dim);
    for(int r = 0; r < 10000; r++) {

      std::vector<std::vector<uint8_t> > datasets;
      for(size_t i=0; i<dim; i++) {
	std::vector<uint8_t> row;
	for(size_t j=0; j<dim; j++) {
	  row.push_back(rng->roll()._poly);
	}
	datasets.push_back(row);
      }
     
    std::vector<size_t> classes;
	classes.push_back((size_t) dim);
  
      nc_generator g(rng.get(), &datasets, classes);
      nc_absorber decoder(rng.get(), dim);

      bool fin = false;
      size_t expected_rank = 0;
      for(size_t i=0; !fin; i++) {
	//      REQUIRE ( decoder.get_current_degree() == expected_rank );
	fin = decoder.inject(g.new_zerotail_linear_combination(dim-(i < dim ? i : 0)));

	//      ++expected_rank;
	//      REQUIRE ( decoder.get_current_degree() == expected_rank );
	//      REQUIRE ( fin == (expected_rank == dim) );
      }

      auto result = decoder.solve();

      CAPTURE(dim);
      CAPTURE(r);
      // check result
      for(size_t i=0; i<datasets.size(); i++) {
	for(size_t j=0; j<datasets.at(0).size(); j++) {
	  REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
	}
      }
    }
  }
}


TEST_CASE("Tests various partial decodings and crazy counter elimination, forward zerotail progression with randomized redundant tails in between", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  for(size_t dim = 5; dim < 20; dim++) {
    for(size_t r = 0; r < 100; r++) {
    
	  std::vector<std::vector<uint8_t> > datasets;
      for(size_t i=0; i<dim; i++) {
		std::vector<uint8_t> row;
		for(size_t j=0; j<dim; j++) {
		  row.push_back(rng->roll()._poly);
		}
		datasets.push_back(row);
      }
	  std::vector<size_t> classes;
	  classes.push_back((size_t) dim);
  
      nc_generator g(rng.get(), &datasets, classes);
      nc_absorber decoder(rng.get(), dim);

      bool fin;
      size_t expected_rank = 0;
      for(size_t i=0; i<dim; i++) {
	REQUIRE ( decoder.get_current_degree() == expected_rank );
	fin = decoder.inject(g.new_zerotail_linear_combination(i+1));

	++expected_rank;
	REQUIRE ( decoder.get_current_degree() == expected_rank );
	// REQUIRE ( decoder.get_maximum_solvable() == expected_rank );
	REQUIRE ( fin == (expected_rank == dim) );

	// inject redundant lower-or-equal ranked zerotail linear combination with 25% chance
	if(rng->roll()._poly < 256/4) {
	  uint8_t zt = rng->roll()._poly % (i+1);
	  fin = decoder.inject(g.new_zerotail_linear_combination(zt+1));
	}
      
      }

      auto result = decoder.solve();

      // check result
      for(size_t i=0; i<datasets.size(); i++) {
	for(size_t j=0; j<datasets.at(0).size(); j++) {
	  REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
	}
      }
    }
  }
}

TEST_CASE("Tests various partial decodings and counter elimination, uses completely random zerotail ranks", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  size_t runs = 100;
  size_t dims = 20;
  
  for(size_t dim = 2; dim < dims; dim++) {
    for(size_t r = 0; r < runs; r++) {
      // printf("checking dimension %zu/%zu, run %zu/%zu:\n", dim+1, dims, r+1, runs);

      std::vector<std::vector<uint8_t> > datasets;
      for(size_t i=0; i<dim; i++) {
		std::vector<uint8_t> row;
		for(size_t j=0; j<dim; j++) {
		row.push_back(rng->roll()._poly);
		}
		datasets.push_back(row);
      }

	  std::vector<size_t> classes;
	  classes.push_back((size_t) dim);
  
      nc_generator g(rng.get(), &datasets, classes);
      nc_absorber decoder(rng.get(), dim);

      bool fin;
      for(size_t i=0;; i++) {
	uint8_t ztr = (rng->roll()._poly % dim) + 1;
	// printf("inserting ztr %u/%zu:", ztr, dim);
	uint16_t old_deg = decoder.get_current_degree();
	fin = decoder.inject(g.new_zerotail_linear_combination(ztr));

	uint16_t new_deg = decoder.get_current_degree();
	if(new_deg > old_deg) {
	  // printf("...innovative (new degree %u)\n", new_deg);
	}
	else {
	  // printf("...non innovative\n");
	}
	if(fin) {
	  break;
	}
      }

      //      printf("solving\n");
      auto result = decoder.solve();

      // check result
      CAPTURE( dim );
      CAPTURE( r );
      for(size_t i = 0; i < dim; i++) {
	for(size_t j = 0; j < dim; j++) {
	  REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
	}
      }
    }
  }
}

TEST_CASE("Tests various partial decodings using naive approach", "[nc_naive]")
{
  auto rng = std::make_shared<simple_rng>();

  size_t runs = 100;
  size_t dims = 20;
  
  for(size_t dim = 2; dim < dims; dim++) {
    for(size_t r = 0; r < runs; r++) {
      // printf("checking dimension %zu/%zu, run %zu/%zu:\n", dim+1, dims, r+1, runs);

      std::vector<std::vector<uint8_t> > datasets;
      for(size_t i=0; i<dim; i++) {
		std::vector<uint8_t> row;
		for(size_t j=0; j<dim; j++) {
		  row.push_back(rng->roll()._poly);
		}
		datasets.push_back(row);
      }

	  std::vector<size_t> classes;
	  classes.push_back((size_t) dim);
  
      nc_generator g(rng.get(), &datasets, classes);
      std::vector<size_t> v(dim-1);
      std::iota (std::begin(v), std::end(v), 2);
      naive_nc_absorber decoder(rng.get(), v);

      bool fin;
      for(size_t i=0;; i++) {
	uint8_t ztr = (rng->roll()._poly % dim) + 1;
	// printf("inserting ztr %u/%zu:", ztr, dim);
	uint16_t old_deg = decoder.get_current_degree();
	fin = decoder.inject(g.new_zerotail_linear_combination(ztr));

	uint16_t new_deg = decoder.get_current_degree();
	if(new_deg > old_deg) {
	  // printf("...innovative (new degree %u)\n", new_deg);
	}
	else {
	  // printf("...non innovative\n");
	}
	if(fin) {
	  break;
	}
      }

      //      printf("solving\n");
      // auto result = decoder.solve();

      // // check result
      // CAPTURE( dim );
      // CAPTURE( r );
      // for(size_t i = 0; i < dim; i++) {
      // 	for(size_t j = 0; j < dim; j++) {
      // 	  REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
      // 	}
      // }
    }
  }
}

TEST_CASE("Reproduces a specific error with 3x3 matrices where two ztr 2 insertions happen", "[nc_partial]")
{
  auto rng = std::make_shared<simple_rng>();

  size_t dim = 3;
  for(size_t r = 0; r < 100; r++) {

      std::vector<std::vector<uint8_t> > datasets;
      for(size_t i=0; i<dim; i++) {
		std::vector<uint8_t> row;
		for(size_t j=0; j<dim; j++) {
		  row.push_back(rng->roll()._poly);
		}
		datasets.push_back(row);
      }

	  std::vector<size_t> classes;
	  classes.push_back((size_t) dim);
  
      nc_generator g(rng.get(), &datasets, classes);
      nc_absorber decoder(rng.get(), dim);

      bool fin;
      for(uint8_t ztr : {2,2,1,3}) {
	// printf("inserting ztr %u/%zu:", ztr, dim);
	uint16_t old_deg = decoder.get_current_degree();
	fin = decoder.inject(g.new_zerotail_linear_combination(ztr));

	uint16_t new_deg = decoder.get_current_degree();
	if(new_deg > old_deg) {
	  // printf("...innovative (new degree %u)\n", new_deg);
	}
	else { 
	  // printf("...non innovative\n");
	}
	if(fin) {
	  break;
	}
      }

      auto result = decoder.solve();

      // check result
      for(size_t i = 0; i < dim; i++) {
	for(size_t j = 0; j < dim; j++) {
	  REQUIRE( (int) datasets.at(i).at(j) == (int) result.at(i).at(j) );
	}
      }
  }
}

TEST_CASE("naive_absorber", "[nc_class_selection_test1]")
{
  simple_rng rng;

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  datasets.push_back({99,12,231});
  datasets.push_back({18,13,17});
  datasets.push_back({32,94,88});
  datasets.push_back({251,5,3});
  datasets.push_back({27,16,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 2);
  classes.push_back((size_t) 4);
  classes.push_back((size_t) 7);
  classes.push_back((size_t) 10);
  
  nc_generator g(&rng, &datasets, classes);

  std::vector<linear_combination> cs;
  for(size_t i=0; i<10; i++) {
    auto c = g.new_zerotail_linear_combination(10);
    cs.push_back(c);
  }

  // decode!
  naive_nc_absorber decoder(&rng, classes );
  for(size_t i=0; i<10; i++) {
    bool fin = decoder.inject(cs.at(i));
    if(i < 9) {
      REQUIRE ( fin == false );
    }
    else {
      REQUIRE ( fin == true );
    }
    // printf("i,fin: %u,%u\n",i,fin);
  }
  auto result = decoder.get_current_degrees();

  for(size_t i=0; i<classes.size(); i++) {
    //REQUIRE( classes.at(i) == result.at(i) );
  }
  
}

TEST_CASE("Classchoice algorithm generator", "[nc_class_selection_test2]")
{
  simple_rng rng;

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  datasets.push_back({99,12,231});
  datasets.push_back({18,13,17});
  datasets.push_back({32,94,88});
  datasets.push_back({251,5,3});
  datasets.push_back({27,16,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 2);
  classes.push_back((size_t) 4);
  classes.push_back((size_t) 7);
  classes.push_back((size_t) 10);
  
  nc_generator g(&rng, &datasets, classes);
  REQUIRE( g.classChoice(4) == 0);
  
  std::vector<uint16_t> fb1;
  fb1.push_back(2);
  fb1.push_back(4);
  fb1.push_back(7);
  fb1.push_back(7);
  g.saveFB((uint32_t)1, fb1);
  REQUIRE( g.classChoice(4) == 3); 

  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=0;
  
  g.saveFB((uint32_t)2, fb1);
  REQUIRE( g.classChoice(4) == 0);
  fb1.at(0)=2;
  fb1.at(1)=2;
  fb1.at(2)=2;
  fb1.at(3)=2;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb1);
  REQUIRE(g.classChoice(4) ==1);
  
  //case1
  std::vector<uint16_t> fb2;
  fb2.push_back(2);
  fb2.push_back(3);
  fb2.push_back(4);
  fb2.push_back(9);
  
  fb1.at(0)=1;
  fb1.at(1)=3;
  fb1.at(2)=4;
  fb1.at(3)=4;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  REQUIRE(g.classChoice(4)== 1);
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  REQUIRE(g.classChoice(1)==0);
  REQUIRE(g.classChoice(1)==2);
  REQUIRE(g.classChoice(3)==2);
  REQUIRE(g.classChoice(3)==3);
  g.delFB((uint32_t)1);
  g.delFB((uint32_t)2);
  REQUIRE(g.classChoice(3)==0);
   
  //
  std::vector<uint16_t> fb3;
  fb3.push_back(0);
  fb3.push_back(0);
  fb3.push_back(7);
  fb3.push_back(7);

  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=9;
  fb2.at(0)=0;
  fb2.at(1)=0;
  fb2.at(2)=7;
  fb2.at(3)=9;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  g.saveFB((uint32_t)3,fb3);
  //send layer 3 instead of layer 2 (count from 0-3)
  REQUIRE(g.classChoice(4)==3);
  
  g.checkComplete();
  REQUIRE(g.completed == false);  
  
  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=10;
  fb2.at(0)=0;
  fb2.at(1)=0;
  fb2.at(2)=7;
  fb2.at(3)=10;
  fb3.at(0)=2;
  fb3.at(1)=4;
  fb3.at(2)=7;
  fb3.at(3)=10;    
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  g.saveFB((uint32_t)3,fb3);
  
  g.checkComplete();
  REQUIRE(g.completed == true);      
}

TEST_CASE("absorber degrees, inject, new combination", "[nc_class_selection_test3]")
{
  simple_rng rng;

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  datasets.push_back({99,12,231});
  datasets.push_back({18,13,17});
  datasets.push_back({32,94,88});
  datasets.push_back({251,5,3});
  datasets.push_back({27,16,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 2);
  classes.push_back((size_t) 4);
  classes.push_back((size_t) 7);
  classes.push_back((size_t) 10);
  
  //new absorber
  naive_nc_absorber decoder(&rng, classes);
  auto feedback= decoder.get_current_degrees();
  for(int i =0; i< (int) feedback.size(); i++){
	  REQUIRE(feedback.at(i) ==0);
  }
  nc_generator g(&rng, &datasets, classes);
  auto lc= g.new_zerotail_linear_combination(7);
  decoder.inject(lc);
  
  feedback= decoder.get_current_degrees();
  for(int i =0; i< (int) feedback.size(); i++){
	  if(i < 2){
		REQUIRE(feedback.at(i) ==0);
	 }
	 else{
		REQUIRE(feedback.at(i) ==1);
	 }
  }
  lc= decoder.new_zerotail_linear_combination(0, 10, 3);
  decoder.inject(lc);
  feedback= decoder.get_current_degrees();
  for(int i =0; i< (int) feedback.size(); i++){
	  if(i < 2){
		REQUIRE(feedback.at(i) ==0);
	 }
	 else{
		REQUIRE(feedback.at(i) ==1);
	 }
  }
}


TEST_CASE("Classchoice algorithm absorber", "[nc_class_selection_test3]")
{
  simple_rng rng;

  std::vector<std::vector<uint8_t> > datasets;
  datasets.push_back({99,123,231});
  datasets.push_back({18,12,17});
  datasets.push_back({32,99,88});
  datasets.push_back({251,9,3});
  datasets.push_back({27,88,65});
  datasets.push_back({99,12,231});
  datasets.push_back({18,13,17});
  datasets.push_back({32,94,88});
  datasets.push_back({251,5,3});
  datasets.push_back({27,16,65});
  
  std::vector<size_t> classes;
  classes.push_back((size_t) 2);
  classes.push_back((size_t) 4);
  classes.push_back((size_t) 7);
  classes.push_back((size_t) 10);
  
  nc_generator generator(&rng, &datasets, classes);
  std::vector<linear_combination> cs;
  for(size_t i=0; i<10; i++) {
    auto c = generator.new_linear_combination();
    cs.push_back(c);
  }
  
  naive_nc_absorber g(&rng, classes);
  for(size_t i=0; i<10; i++) {
    bool fin = g.inject(cs.at(i));
  }

  REQUIRE( g.classChoice(4) == 0);
  
  std::vector<uint16_t> fb1;
  fb1.push_back(2);
  fb1.push_back(4);
  fb1.push_back(7);
  fb1.push_back(7);
  g.saveFB((uint32_t)1, fb1);
  REQUIRE( g.classChoice(4) == 3); 

  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=0;
  
  g.saveFB((uint32_t)2, fb1);
  REQUIRE( g.classChoice(4) == 0);
  fb1.at(0)=2;
  fb1.at(1)=2;
  fb1.at(2)=2;
  fb1.at(3)=2;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb1);
  REQUIRE(g.classChoice(4) ==1);
  
  //case1
  std::vector<uint16_t> fb2;
  fb2.push_back(2);
  fb2.push_back(3);
  fb2.push_back(4);
  fb2.push_back(9);
  
  fb1.at(0)=1;
  fb1.at(1)=3;
  fb1.at(2)=4;
  fb1.at(3)=4;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  REQUIRE(g.classChoice(4)== 1);
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  REQUIRE(g.classChoice(1)==0);
  REQUIRE(g.classChoice(1)==2);
  REQUIRE(g.classChoice(3)==2);
  REQUIRE(g.classChoice(3)==3);
  g.delFB((uint32_t)1);
  g.delFB((uint32_t)2);
  REQUIRE(g.classChoice(3)==0);

std::vector<uint16_t> fb3;
  fb3.push_back(0);
  fb3.push_back(0);
  fb3.push_back(7);
  fb3.push_back(7);

  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=9;
  fb2.at(0)=0;
  fb2.at(1)=0;
  fb2.at(2)=7;
  fb2.at(3)=9;
  
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  g.saveFB((uint32_t)3,fb3);
  //send layer 3 instead of layer 2 (count from 0-3)
  REQUIRE(g.classChoice(4)==3);
  
  g.checkComplete();
  REQUIRE(g.completed == false);  
  
  
  fb1.at(0)=0;
  fb1.at(1)=0;
  fb1.at(2)=0;
  fb1.at(3)=10;
  fb2.at(0)=0;
  fb2.at(1)=0;
  fb2.at(2)=7;
  fb2.at(3)=10;
  fb3.at(0)=2;
  fb3.at(1)=4;
  fb3.at(2)=7;
  fb3.at(3)=10;    
  g.saveFB((uint32_t)1,fb1);
  g.saveFB((uint32_t)2,fb2);
  g.saveFB((uint32_t)3,fb3);
  
  g.checkComplete();
  REQUIRE(g.completed == true);    
}
