def reliabilityContinuo(new2, top, bottom):
    'take in input the output of calculation by vera continuo'

    reliability = []
    for i in new2.index:
        if isinstance(new2['Grouping Value'][i], float) and isinstance(new2['Knn'][i], float):
            delta_GK = round(abs(new2['Grouping Value'][i]-new2['Knn'][i])*100/abs(top-bottom), 3)
        else: delta_GK = 100
        
        if isinstance(new2['Knn'][i], float) and isinstance(new2['LocalModel'][i], float):
            delta_KL = round(abs(new2['Knn'][i]-new2['LocalModel'][i])*100/abs(top-bottom), 3)
        else: delta_KL = 100
            
        if isinstance(new2['Grouping Value'][i], float) and isinstance(new2['LocalModel'][i], float):
            delta_GL = round(abs(new2['Grouping Value'][i]-new2['LocalModel'][i])*100/abs(top-bottom), 3)
        else: delta_GL = 100

        if isinstance(new2['err'][i], str):
            new2['err'][i] = 100
        
        if new2['reasoning'][i] == 'no': # i termini per il knn sono buoni
            
            if new2['LocalModel'][i] == 'no': # il grouping è buono 
                if delta_GK<20:
                    reliability.append('Good reliability') 
                elif delta_KL<20 and new2['err'][i]<0.51: 
                    reliability.append('Good reliability')
                else: reliability.append('Moderate reliability')
            else: # i termini per il grouping non sono buoni          
                if delta_KL<15 and new2['err'][i]<=0.50:
                    reliability.append('Good reliability')
                else:
                    reliability.append('Moderate reliability')                    
        else: # termini knn non buoni
            if new2['LocalModel'][i] == 'no': # termini del grouping sono buoni
                if delta_GL<20 and new2['err'][i]<=0.50: reliability.append('Moderate reliability')
                else: reliability.append('Low reliability')
            else: # termini del grouping sono buoni e del knn
                if  new2['err'][i]<0.20: reliability.append('Moderate reliability')
                else: reliability.append('Low reliability')
    
    return reliability