#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

std::string getRandomThankYou() {
    std::vector<std::string> thankYouMessages = {
        "Thank you for your support!",
        "You're awesome! Enjoy your beer!",
        "Cheers to you!",
        "Thanks for the coffee!",
        "You just made my day!",
        "HOPL appreciates your kindness!",
        "Your support means a lot!",
        "Enjoy your drink, and thanks again!"
    };

    int index = rand() % thankYouMessages.size();
    return thankYouMessages[index];
}

int main() {
    srand(static_cast<unsigned int>(time(0)));

    std::cout << "Welcome to the HOPL License Project!" << std::endl;
    std::cout << "If you enjoyed using this software, consider supporting the author with a beer or coffee!" << std::endl;
    std::cout << "Press Enter to receive a random thank you message..." << std::endl;

    std::cin.get();

    std::string message = getRandomThankYou();
    std::cout << "\n" << message << "\n" << std::endl;

    return 0;
}
