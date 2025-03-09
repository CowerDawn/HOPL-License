#!/bin/bash

echo "Welcome to the HOPL Donation Script!"

donation_amount=0

while true; do
    echo -n "Would you like to donate to the author? (y/n): "
    read choice

    if [[ $choice == "y" ]]; then
        echo -n "Enter donation amount (in USD): "
        read amount
        donation_amount=$((donation_amount + amount))
        echo "Thank you for your donation of $amount USD! Total donations: $donation_amount USD"
    elif [[ $choice == "n" ]]; then
        echo "Thank you for your support! Total donations: $donation_amount USD"
        break
    else
        echo "Invalid choice. Please enter 'y' or 'n'."
    fi
done

