import sqlite3
from stellarutil.simulation import Simulation

# Get list of stars
m10r = Simulation('m10r_res250md', species=['star'])

# Connect to the Database
conn = sqlite3.connect('star_data.db')  # Replace with your database name
cursor = conn.cursor()


# Create a table if it doesn't exist
cursor.execute('''
    CREATE TABLE IF NOT EXISTS stars (
        id INTEGER PRIMARY KEY,
        x REAL,
        y REAL,
        z REAL,
        vx REAL,
        vy REAL,
        vz REAL,
        a REAL,
        r REAL
    )
''')



# Iterate Through the Array and Insert Data
for index in range(100):
    for star in m10r.get_halo(index).stars:
        cursor.execute('''
            INSERT INTO stars (x, y, z, vx, vy, vz, a, r)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (float(star.x), float(star.y), float(star.z), float(star.vx), float(star.vy), float(star.vz), float(star.a), float( (star.x**2 + star.y**2)**0.5) ))
            

# Commit the changes and close the connection
conn.commit()
conn.close()
