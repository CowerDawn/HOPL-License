let activeUsers = 0;

function updateActiveUsers() {
    document.getElementById('activeUsers').textContent = activeUsers;
}

setInterval(() => {

    activeUsers += Math.floor(Math.random() * 3);
    updateActiveUsers();
}, 5000);

setInterval(() => {
    if (activeUsers > 0) {
        activeUsers -= Math.floor(Math.random() * 2);
        updateActiveUsers();
    }
}, 10000);


updateActiveUsers();
