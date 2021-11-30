# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 10:59:18 2021

@author: Dario Ledesma
"""

import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
class radiationClass:
    def __init__(self, dataFrame, LAT, LONG, ALT):
        self.dataFrame = dataFrame
        self.LAT= LAT
        self.LONG = LONG
        self.ALT = ALT
        self.GMT = 0
        self.dataFrame['Fecha'] = pd.to_datetime(self.dataFrame['Fecha'])
        self.dataFrame['Year'] = self.dataFrame['Fecha'].dt.year
        self.dataFrame['Dia juliano'] = self.dataFrame['Fecha'].dt.day_of_year
        self.dataFrame['Mes'] = self.dataFrame['Fecha'].dt.month
        self.dataFrame['Hora'] = self.dataFrame['Fecha'].dt.time
        self.dataFrame['Hora reloj'] = self.dataFrame['Fecha'].dt.hour + self.dataFrame['Fecha'].dt.minute * 0.0166667
        self.dataFrame['Angulo diario'] = self.daily_angle(self.dataFrame['Dia juliano'])
        self.dataFrame['Ecuacion del tiempo'] = self.generate_time_equation(self.dataFrame['Angulo diario'])
        self.dataFrame['Hora solar'] = self.generate_hora_solar(self.dataFrame['Hora reloj'], self.dataFrame['Ecuacion del tiempo']  )
        self.dataFrame['Declinacion'] = self.generate_declination(self.dataFrame['Dia juliano'])
        self.dataFrame['Declimacion en radianes'] = np.multiply(self.dataFrame["Declinacion"],  0.017453)
        self.dataFrame['Angulo horario'] = self.generate_hour_angle(self.dataFrame['Hora solar'])
        self.dataFrame['Angulo horario en radianes'] = np.multiply(self.dataFrame["Angulo horario"],  0.017453)
        self.dataFrame['E0'] = self.generate_e0(self.dataFrame['Dia juliano'])
        self.dataFrame['Cos tita z'] = self.generate_cos_tita_z(self.dataFrame['Declimacion en radianes'], self.dataFrame["Angulo horario en radianes"])
        self.dataFrame['Tita z'] = self.generate_tita_z()
        self.dataFrame['Altura solar'] = 90 - self.dataFrame['Tita z'] 
        self.dataFrame['Irr TOA wm2'] = self.generate_irradiancia_ext(self.dataFrame['Cos tita z'] , 1367.00, self.dataFrame['E0'])
        self.ktr = 0.7 + 1.6391 * math.pow(10,-3) * math.pow(self.ALT, 0.5500)
        self.dataFrame['Cos tita z grad'] = self.dataFrame['Cos tita z'].apply(math.degrees)
        self.dataFrame['MA'] = self.generate_ma()
        self.dataFrame['GHIcc'] = self.generate_GHIcc()
        # self.grupoDiario = self.day_mean()
        # self.grupoYear = self.year_mean()
        
    def daily_angle(self, day):
        return np.divide(np.multiply(day-1, 2.0 * math.pi),365)
    



    def generate_time_equation(self, angles):
        return  (0.0000075+0.001868* angles.apply(math.cos)-0.032077* angles.apply(math.sin)-0.014615* (2*angles).apply(math.cos) -0.04089*  (2*angles).apply(math.sin))*229.2



    def generate_hora_solar(self, clock_hour, ecuation_time):
        A = 1;
        if(self.GMT<=0):
            A = -1;
        return clock_hour + (4*((A*15 * self.GMT)-(A*self.LONG))+ ecuation_time)/60

    def generate_declination(self, julian_days):
        return 23.45 * ((360*(284 + julian_days)).apply(math.radians)/365).apply(math.sin)



    def generate_hour_angle(self, solar_hour):
        new_values = np.multiply(12-solar_hour, 15)
        return new_values
    
    def generate_cos_tita_z(self, declination_rad, hour_angle_rad ):  
        C2 = 0.017453 * self.LAT
        L2 = declination_rad
        N2 = hour_angle_rad
        return (math.cos(C2)* L2.apply(math.cos) * N2.apply(math.cos))+(math.sin(C2)* L2.apply(math.sin))
    
    
    def generate_e0(self, days):
        return (1+0.033* (2* math.pi * days /365).apply(math.cos) )
    
    def generate_irradiancia_ext(self, cos_tita_z, TSI, E0):
        return np.where(cos_tita_z<0,0, TSI * E0 * cos_tita_z )
    
    def generate_tita_z(self):
        return self.dataFrame['Cos tita z'].apply(math.acos).apply(math.degrees)
    
    
    def generate_ma(self):
        
        med = 0.15 *  np.power(93.885 - self.dataFrame['Cos tita z grad'], -1.253)
        return 1 / (self.dataFrame['Cos tita z grad'] + med)
        
        #1 / (self.dataFrame['Cos tita z grad'] + 0.15* math.pow(93.885 - self.dataFrame['Cos tita z grad'], -1.253))
    
    
    def generate_GHIcc(self):
        med = np.power(self.dataFrame['MA'], 0.678)
        return self.dataFrame['Irr TOA wm2'] * np.power(self.ktr, med)
    
    
    def export(self, name):
        self.dataFrame.to_csv(name, sep=",")
    
    def day_mean(self):
        data = self.dataFrame[['Dia juliano', 'Clear sky GHI', 'Clear sky BHI', 'Clear sky DHI', 'Clear sky DNI', 'Hora reloj', 'Year']]
        return data.groupby(['Hora reloj']).mean()
    
    
    def year_mean(self):
        data = self.dataFrame[['Dia juliano', 'Clear sky GHI', 'Clear sky BHI', 'Clear sky DHI', 'Clear sky DNI', 'Hora reloj', 'Year']]
        return data.groupby(['Dia juliano', 'Hora reloj']).mean().reset_index()
    
    
    def month_mean(self):
         data = self.dataFrame[['Clear sky GHI', 'Clear sky DNI','Clear sky DHI','Hora reloj', 'Dia juliano' ,'Mes','Year']]
         return data.groupby(['Year', 'Hora reloj']).mean().reset_index()
    
    def all_mean(self):
         data = self.dataFrame[['Clear sky GHI', 'Clear sky DNI','Clear sky DHI','Hora reloj', 'Dia juliano' ,'Mes','Year']]
         return data.groupby(['Hora reloj']).mean().reset_index()
    
    def plot_day(self):
        self.grupoDiario = self.grupoDiario.reset_index()
        hora = self.grupoDiario['Hora reloj']
        GHI = self.grupoDiario['Clear sky GHI']
        BHI = self.grupoDiario['Clear sky BHI']
        DHI = self.grupoDiario['Clear sky DHI']
        DNI = self.grupoDiario['Clear sky DNI']
        plt.plot(hora, GHI, label = "GHI")
        plt.plot(hora, BHI, label = "BHI")
        plt.plot(hora, DHI, label = "DHI")
        plt.plot(hora, DNI, label = "DNI")
        plt.legend()
        plt.show()
        
        
    def plot_day_of_year(self, day):
        data = self.dataFrame[self.dataFrame["Dia juliano"] == day]
        hora = data['Hora reloj']
        GHI = data['Global']
        TOA = data['Irr TOA wm2']
        GHIcc = data['GHIcc']
        plt.plot(hora, GHI, label = "GHI")
        plt.plot(hora, TOA, label = "TOA")
        plt.plot(hora, GHIcc, label = "GHIcc")
        plt.legend()
        plt.ylim([0, 2000])
        plt.title("dia "+ str(day))
        plt.show()
        
        
    def plot_solar_map(self, day):
        data = self.dataFrame[self.dataFrame["Dia juliano"] == day]
        hora = data['Tita z']
        altura = data['Altura solar']
        
        plt.plot(hora, altura, label = "Diagrama solar")
        plt.legend()
        #plt.ylim([0, 2000])
        plt.title("dia "+ str(day))
        plt.show()
    
    
    def plot_all_years(self):
        fig, (axs1, axs2) = plt.subplots(1, 2)
        #fig.subplots_adjust(hspace=0.4, wspace=0.4)
        #ax = fig.add_subplot(113)
        #ax[0,0].text(0.5, 0.5, str('Año, 2007' ),fontsize=10, ha='center')
        data = self.grupoMes[self.grupoMes["Year"] == 2007]
        dia = data['Hora reloj']
        GHI = data['Clear sky GHI']
        BHI = data['Clear sky BHI']
        axs1.plot(dia,GHI, label="GHI")
        axs1.plot(dia,BHI, label="BHI")
        
        
        #ax = fig.add_subplot(111)
        #ax[0,1].text(0.5, 0.5, str('Año, 2008' ),fontsize=10, ha='center')
        data = self.grupoMes[self.grupoMes["Year"] == 2020]
        dia = data['Hora reloj']
        GHI = data['Clear sky GHI']
        BHI = data['Clear sky BHI']
        axs2.plot(dia,GHI, label="GHI")
        axs2.plot(dia,BHI, label="BHI")
        
        
        plt.show()
    
    def plot_year(self):
        years = self.grupoMes['Year'].unique()
        # fig, axs = plt.subplots(2, int(len(years)/2)
        # fig.suptitle('Vertically stacked subplots')
        
        fig = plt.figure()
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        for i in range(1, 20):
            
            
            
            ax = fig.add_subplot(4, 5, i)
            #ax.text(0.5, 0.5, str((4, 5, i)),fontsize=18, ha='center')
            #ax.text(0.5, 0.5, str(("Año", )),fontsize=18, ha='center')
            
            if(i<=len(years)):
                
                year = years[i-1]
                ax.text(0.5, 0.5, str(("Año", year )),fontsize=18, ha='center')
                data = self.grupoMes[self.grupoMes["Year"] == year]
                dia = data['Hora reloj']
                GHI = data['Clear sky GHI']
                BHI = data['Clear sky BHI']
                ax.plot(dia,GHI, label="GHI")
                ax.plot(dia,BHI, label="BHI")
            plt.show()
        
        # for index , val in enumerate(years):
        #     data = self.grupoMes[self.grupoMes["Year"] == val]
        #     dia = data['Hora reloj']
        #     GHI = data['Clear sky GHI']
        #     BHI = data['Clear sky BHI']
        #     plt.plot(dia, GHI, label = "GHI")
        #     axs[index,1].plot(dia, GHI, label = "GHI")
        #     axs[index,0].plot(dia, BHI, label = "BHI")



# yuto = pd.read_csv('EEA_Yuto_31082021_1236-489.csv', sep=",", usecols=[2,3,7,11])
# yuto.columns = ['Fecha', 'Global', 'Directa', 'Difusa']

# yuto['Fecha'] = pd.to_datetime(yuto['Fecha'])
# yuto = yuto.sort_values(by='Fecha')

# yuto.to_csv('yuto.csv', index=False)


# abraPampa = abraPampa.set_index(abraPampa['Fecha'])

# abraPampa.index = pd.DatetimeIndex(abraPampa.index)
# abraPampa = abraPampa.reindex(r, fill_value='NaN')


# abraPampa['Fecha'] = abraPampa.index


# abraPampa = pd.read_csv('abraPampa_2018.csv', sep=",")
# abraPampa = radiationClass(abraPampa, LAT=-65.69, LONG=-22.72)

#abraPampa = radiationClass(abraPampa, LAT=-65.69, LONG=-65.69 ).dataFrame



# yuto = pd.read_csv('yuto.csv', sep=",")
# yuto['Fecha'] = pd.to_datetime(yuto['Fecha'])
# yuto2018 = yuto[yuto['Fecha'].dt.year == 2018]
# r = pd.date_range(start='2018-01-01', end = '2019-01-01', freq='min')
# r = r[:-1]

# yuto2018= yuto2018.set_index(yuto2018['Fecha'])

# yuto2018.index = pd.DatetimeIndex(yuto2018.index)
# yuto2018 = yuto2018.reindex(r, fill_value='NaN')

# yuto2018['Fecha'] = yuto2018.index
# yuto2018.to_csv('yuto_2018.csv', index=False)



salta = pd.read_csv('salta_2018.csv', sep=",")
abrapampa = pd.read_csv('abraPampa_2018.csv', sep=",")
yuto = pd.read_csv('yuto_2018.csv', sep=",")


salta_geo = radiationClass(salta, -24.78, -65.41, 1152)
abrapampa_geo = radiationClass(abrapampa, -22.72, -65.69,3487)
yuto_geo = radiationClass(yuto, -23.38, -64.28, 340)

salta_geo.plot_day_of_year(6)



