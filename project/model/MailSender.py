import smtplib

class MailSender:


    def sendMail(self):
        remitente = "Haplotype Searcher <paulavillanueva.93@gmail.com>"
        destinatario = "Paula Villanueva <paulavillanueva.22@gmail.com>"
        asunto = "E-mal HTML enviado desde Python"
        mensaje = """Hola!<br/> <br/> 
        Este es un <b>e-mail</b> enviando desde <b>Python</b> 
        """

        email = """From: %s 
        To: %s 
        MIME-Version: 1.0 
        Content-type: text/html 
        Subject: %s 
        
        %s
        """ % (remitente, destinatario, asunto, mensaje)
        try:
            server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
            server.ehlo()
            server.login('paulavillanueva.93@gmail.com', 'panaholma')
            server.sendmail(remitente, destinatario, email)
            print("Correo enviado")
        except Exception as e :
            print(e)
