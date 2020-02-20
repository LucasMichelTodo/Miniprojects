import mysql.connector as mysql

db = mysql.connect(
    host="localhost",
    user="root",
    passwd="lucas888"
)

cursor = db.cursor()

cursor.execute("CREATE DATABASE datacamp")
cursor.execute("SHOW DATABASES")
databases = cursor.fetchall()
for d in databases:
    print(d)

db = mysql.connect(host="localhost",
                   user="root",
                   passwd="lucas888",
                   database="datacamp")

cursor = db.cursor()
cursor.execute(
    "CREATE TABLE users (name VARCHAR(255), user_name VARCHAR(255))")

# getting all the tables which are present in 'datacamp' database
cursor.execute("SHOW TABLES")

tables = cursor.fetchall()  # it returns list of tables present in the database

# showing all the tables one by one
for table in tables:
    print(table)

cursor.execute("DROP TABLE users")

cursor.execute(
    "CREATE TABLE users (id INT(11) NOT NULL AUTO_INCREMENT PRIMARY KEY, name VARCHAR(255), user_name VARCHAR(255))")

# 'DESC table_name' is used to get all columns information
cursor.execute("DESC users")

# it will print all the columns as 'tuples' in a list
print(cursor.fetchall())

query = "INSERT INTO users (name, user_name) VALUES (%s, %s)"
values = ("Lucas", "Gnaaaarg")

cursor.execute(query, values)
db.commit()
print(cursor.rowcount, "record inserted")

values = [("Peter", "peter"),
          ("Amy", "amy"),
          ("Michael", "michael"),
          ("Hennah", "hennah")]

cursor.executemany(query, values)
db.commit()
print(cursor.rowcount, "record inserted")

cursor.execute("SELECT * FROM users")
records = cursor.fetchall()

for r in records:
    print(r)

query = "SELECT user_name FROM users"
cursor.execute(query)

usernames = cursor.fetchall()

query = "SELECT name, user_name FROM users"

usernames
