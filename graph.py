import matplotlib.pyplot as plt

# Чтение данных из файла
with open("qwe.txt") as file:
    data = file.readlines()
    data = [line.strip().split(" ") for line in data]
    #x = [float(line[0]) for line in data]
    y = [float(line[1])  for line in data]
    for line in data:
        print(line[1])
    

plt.title("Амплитуда умноженная на 10")
x = numbers = list(range(0, 32000))
plt.ylim(-180, 180)
# Построение графика
plt.plot(x, y)
#print(y.size)
# Настройка заголовка и меток осей
plt.title("Фаза")

#ax = plt.gca()
##ax.spines['left'].set_position('center')     
#ax.spines['bottom'].set_position('center')
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)


#plt.xlabel("Метка оси x")
#plt.ylabel("Метка оси y") 

# Отображение графика
plt.show()