use rand::Rng;
use std::io;

fn main() {
    let mut rng = rand::thread_rng();
    let messages = vec![
        "Thank you for forking the project!",
        "Your fork is appreciated!",
        "Cheers to your contribution!",
        "Thanks for supporting HOPL!",
        "Your fork makes the project better!",
    ];

    println!("Welcome to the HOPL Fork Simulator!");

    loop {
        println!("\nWould you like to create a fork? (y/n): ");
        let mut input = String::new();
        io::stdin().read_line(&mut input).expect("Failed to read input");

        match input.trim().to_lowercase().as_str() {
            "y" => {
                let message_index = rng.gen_range(0..messages.len());
                println!("{}", messages[message_index]);
            }
            "n" => {
                println!("Thank you for using the HOPL Fork Simulator!");
                break;
            }
            _ => println!("Invalid input. Please enter 'y' or 'n'."),
        }
    }
}
use rand::Rng;
use std::io;

fn main() {
    let mut rng = rand::thread_rng();
    let messages = vec![
        "Thank you for forking the project!",
        "Your fork is appreciated!",
        "Cheers to your contribution!",
        "Thanks for supporting HOPL!",
        "Your fork makes the project better!",
    ];

    println!("Welcome to the HOPL Fork Simulator!");

    loop {
        println!("\nWould you like to create a fork? (y/n): ");
        let mut input = String::new();
        io::stdin().read_line(&mut input).expect("Failed to read input");

        match input.trim().to_lowercase().as_str() {
            "y" => {
                let message_index = rng.gen_range(0..messages.len());
                println!("{}", messages[message_index]);
            }
            "n" => {
                println!("Thank you for using the HOPL Fork Simulator!");
                break;
            }
            _ => println!("Invalid input. Please enter 'y' or 'n'."),
        }
    }
}
