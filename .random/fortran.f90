program heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)

    ! Граничные условия (правая граница, например, изолирована)
    u(nx) = u(nx - 1)

    ! Основной цикл по времени
    do j = 1, nt
        ! Вычисление новой температуры
        do i = 2, nx - 1
            u_new(i) = u(i) + alpha * dt / dx**2 * (u(i+1) - 2.0 * u(i) + u(i-1))
        end do

        ! Обновление граничных условий
        u_new(1) = 100.0  ! Левая граница
        u_new(nx) = u_new(nx - 1)  ! Правая граница

        ! Обновление массива температуры
        u = u_new
    end do

    ! Вывод результатов в файл
    open(unit=10, file='temperature_profile.txt', status='replace')
    do i = 1, nx
        write(10, *) (i-1)*dx, u(i)
    end do
    close(10)

    print *, "Результаты сохранены в файл temperature_profile.txt"
end program heat_equationprogram heat_equation
    implicit none
    integer, parameter :: nx = 100       ! Количество точек по пространству
    integer, parameter :: nt = 1000      ! Количество шагов по времени
    real(8), parameter :: alpha = 0.01   ! Коэффициент температуропроводности
    real(8), parameter :: L = 1.0        ! Длина стержня
    real(8), parameter :: dx = L / (nx - 1)  ! Шаг по пространству
    real(8), parameter :: dt = 0.0001    ! Шаг по времени
    real(8), parameter :: T = nt * dt    ! Общее время моделирования
    real(8) :: u(nx), u_new(nx)          ! Массивы для температуры
    integer :: i, j

    ! Начальные условия (например, температура в начальный момент времени)
    u = 0.0
    u(1) = 100.0  ! Левая граница (постоянная температура)
