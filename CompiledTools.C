// C++ includes...
#include <iostream>

// Function to write out a sweet status bar!
void StatusBar(Int_t current, Int_t total, Int_t period){
	if(current != 0){
		for(Int_t j = 0; j < ((total/period)+2); j++){std::cout << "\b";}
	}
	std::cout << "[";
	for(Int_t j = 0; j < (total/period); j++){
		if(j <= (current/period)){std::cout << "*";}else{std::cout << " ";}
	}
	std::cout << "]" << std::flush;
	
	return;
}
