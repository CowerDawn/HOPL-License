package main

import (
    "fmt"
)

func main() {
    var beerCount int
    var choice string

    fmt.Println("Welcome to the HOPL Beer Counter!")

    for {
        fmt.Print("\nWould you like to buy the author a beer? (y/n): ")
        fmt.Scan(&choice)

        if choice == "y" {
            beerCount++
            fmt.Printf("Thank you for the beer! Total beers: %d\n", beerCount)
        } else if choice == "n" {
            fmt.Printf("\nThank you for your support! Total beers: %d\n", beerCount)
            break
        } else {
            fmt.Println("Invalid choice. Please enter 'y' or 'n'.")
        }
    }
}
