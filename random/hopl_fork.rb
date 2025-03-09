messages = [
    "Thank you for forking the project!",
    "Your fork is appreciated!",
    "Cheers to your contribution!",
    "Thanks for supporting HOPL!",
    "Your fork makes the project better!",
]

puts "Welcome to the HOPL Fork Simulator!"

loop do
    print "\nWould you like to create a fork? (y/n): "
    choice = gets.chomp.downcase

    if choice == "y"
        puts messages.sample
    elsif choice == "n"
        puts "Thank you for using the HOPL Fork Simulator!"
        break
    else
        puts "Invalid input. Please enter 'y' or 'n'."
    end
end
