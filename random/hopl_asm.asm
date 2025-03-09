section .data
    message db 'HOPL', 0xA  
    len equ $ - message      

section .text
    global _start

_start:
    mov eax, 1              
    mov edi, 1              
    lea rsi, [message]      
    mov edx, len            
    syscall

    mov eax, 60             
    xor edi, edi            
    syscall
