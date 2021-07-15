import hashlib
import http.client
import json
import os
import random
import urllib

import openpyxl
import pandas as pd
from openpyxl.styles import Font

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)

def translate(text):
    appid = '20210715000889041'  # 填写你的appid
    secretKey = 'dJ6l4YNxso_x195V65HG'  # 填写你的密钥

    httpClient = None
    myurl = '/api/trans/vip/translate'

    fromLang = 'en'   #原文语种
    toLang = 'zh'   #译文语种
    salt = random.randint(32768, 65536)
    q = text
    sign = appid + q + str(salt) + secretKey
    sign = hashlib.md5(sign.encode()).hexdigest()
    myurl = myurl + '?appid=' + appid + '&q=' + urllib.parse.quote(q) + '&from=' + fromLang + '&to=' + toLang + '&salt=' + str(
    salt) + '&sign=' + sign

    try:
        httpClient = http.client.HTTPConnection('api.fanyi.baidu.com')
        httpClient.request('GET', myurl)

        # response是HTTPResponse对象
        response = httpClient.getresponse()
        result_all = response.read().decode("utf-8")
        result = json.loads(result_all)
        trans_result = result['trans_result'][0]['dst']

    except Exception as e:
        print (e)
    finally:
        if httpClient:
            httpClient.close()
    
    return trans_result

def trans_ref(excel_path, sheet_name):
    lis = []
    translatefile = pdb_path + '%s_Translation.xlsx' % sheet_name
    writer = pd.ExcelWriter(translatefile)
    data = pd.read_excel(excel_path, sheet_name=sheet_name)
    for index, row in data.iterrows():
        translation = translate(str(row['Abstract']))
        print('Translation:', translation)
        lis.append({'DOI' : row['DOI'],'Title' : row['Article Title'],'Journal' : row['Source Title'], 'Abstract' : row['Abstract'],'Translation' : translation})

    df = pd.DataFrame(lis)
    df.to_excel(writer, sheet_name)
    writer.save()

    wb = openpyxl.load_workbook(translatefile)
    font = Font(name='等线',size=12)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(translatefile)

if __name__ == '__main__':
    excel_path = input('Excel Path:')
    sheet_name = input('Sheet Name:')
    trans_ref(excel_path, sheet_name)
