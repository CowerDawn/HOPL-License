#include <stdio.h>
#include <stdlib.h>

int main() {
    int beerCount = 0;
    int coffeeCount = 0;
    char choice;

    printf("Welcome to the HOPL Beer & Coffee Counter!\n");

    while (1) {
        printf("\nWhat would you like to buy for the author?\n");
        printf("1. Beer \n");
        printf("2. Coffee \n");
        printf("3. Exit\n");
        printf("Your choice: ");
        scanf(" %c", &choice);

        if (choice == '1') {
            beerCount++;
            printf("Thank you for the beer! Total beers: %d\n", beerCount);
        } else if (choice == '2') {
            coffeeCount++;
            printf("Thank you for the coffee! Total coffees: %d\n", coffeeCount);
        } else if (choice == '3') {
            printf("\nThank you for your support! Here's the summary:\n");
            printf("Total beers: %d\n", beerCount);
            printf("Total coffees: %d\n", coffeeCount);
            break;
        } else {
            printf("Invalid choice. Please try again.\n");
        }
    }

    return 0;
}
